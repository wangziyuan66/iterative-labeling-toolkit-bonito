"""
Bonito Basecaller
"""

import os
import sys
import numpy as np
from tqdm import tqdm
from time import perf_counter
from functools import partial
from datetime import timedelta
from itertools import islice as take
from itertools import permutations
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from bonito.nn import fuse_bn_
from bonito.aligner import align_map, Aligner
from bonito.reader import Reader, ReadChunk
from bonito.io import CTCWriter, Writer, biofmt
from bonito.mod_util import call_mods, load_mods_model
from bonito.cli.download import Downloader, models, __models_dir__
from bonito.multiprocessing import process_cancel, process_itemmap
from bonito.util import column_to_set, load_symbol, load_model, init, tqdm_environ
from pod5_shuffle import get_reads
from sr_io import srCTCWriter

import torch


def main(args):

    init(args.seed, args.device)

    try:
        reader = Reader(args.reads_directory, args.recursive)
        sys.stderr.write("> reading %s\n" % reader.fmt)
    except FileNotFoundError:
        sys.stderr.write("> error: no suitable files found in %s\n" % args.reads_directory)
        exit(1)

    fmt = biofmt(aligned=args.reference is not None)

    if args.reference and args.reference.endswith(".mmi") and fmt.name == "cram":
        sys.stderr.write("> error: reference cannot be a .mmi when outputting cram\n")
        exit(1)
    elif args.reference and fmt.name == "fastq":
        sys.stderr.write(f"> warning: did you really want {fmt.aligned} {fmt.name}?\n")
    else:
        sys.stderr.write(f"> outputting {fmt.aligned} {fmt.name}\n")

    if args.model_directory in models and not (__models_dir__ / args.model_directory).exists():
        sys.stderr.write("> downloading model\n")
        Downloader(__models_dir__).download(args.model_directory)

    sys.stderr.write(f"> loading model {args.model_directory}\n")
    try:
        model = load_model(
            args.model_directory,
            args.device,
            weights=args.weights if args.weights > 0 else None,
            chunksize=args.chunksize,
            overlap=args.overlap,
            batchsize=args.batchsize,
            quantize=args.quantize,
            use_koi=True,
        )
        model = model.apply(fuse_bn_)
    except FileNotFoundError:
        sys.stderr.write(f"> error: failed to load {args.model_directory}\n")
        sys.stderr.write(f"> available models:\n")
        for model in sorted(models): sys.stderr.write(f" - {model}\n")
        exit(1)

    if args.verbose:
        sys.stderr.write(f"> model basecaller params: {model.config['basecaller']}\n")

    basecall = load_symbol(args.model_directory, "basecall")

    mods_model = None
    if args.modified_base_model is not None or args.modified_bases is not None:
        sys.stderr.write("> loading modified base model\n")
        mods_model = load_mods_model(
            args.modified_bases, args.model_directory, args.modified_base_model,
            device=args.modified_device,
        )
        sys.stderr.write(f"> {mods_model[1]['alphabet_str']}\n")

    if args.reference:
        sys.stderr.write("> loading reference\n")
        aligner = Aligner(args.reference, preset=args.mm2_preset, k=5)
        if not aligner:
            sys.stderr.write("> failed to load/build index\n")
            exit(1)
    else:
        aligner = None

    if args.save_ctc and not args.reference:
        sys.stderr.write("> a reference is needed to output ctc training data\n")
        exit(1)

    if fmt.name != 'fastq':
        groups, num_reads = reader.get_read_groups(
            args.reads_directory, args.model_directory,
            n_proc=8, recursive=args.recursive,
            read_ids=column_to_set(args.read_ids), skip=args.skip,
            cancel=process_cancel()
        )
    else:
        groups = []
        num_reads = None

    reads = get_reads(
        args.reads_directory, n_proc=8, recursive=args.recursive,
        read_ids=column_to_set(args.read_ids), skip=args.skip,
        do_trim=not args.no_trim,
        scaling_strategy=model.config.get("scaling"),
        norm_params=(model.config.get("standardisation")
                     if (model.config.get("scaling") and
                         model.config.get("scaling").get("strategy") == "pa")
                     else model.config.get("normalisation")
                     ),
        cancel=process_cancel()
    )

    if args.verbose:
        sys.stderr.write(f"> read scaling: {model.config.get('scaling')}\n")
    
    if args.max_reads:
        reads = take(reads, args.max_reads)
        if num_reads is not None:
            num_reads = min(num_reads, args.max_reads)

    if args.save_ctc:
        reads = (
            chunk for read in reads
            for chunk in read_chunks(
                read,
                chunksize=model.config["basecaller"]["chunksize"],
                overlap=model.config["basecaller"]["overlap"]
            )
        )
        if args.sr or args.chemo:
            ResultsWriter = srCTCWriter
        else:
            ResultsWriter = CTCWriter
    else:
        ResultsWriter = Writer

    results = basecall(
        model, reads, reverse=args.revcomp, rna=args.rna,
        batchsize=model.config["basecaller"]["batchsize"],
        chunksize=model.config["basecaller"]["chunksize"],
        overlap=model.config["basecaller"]["overlap"]
    )
    
    if mods_model is not None:
        if args.modified_device:
            results = ((k, call_mods(mods_model, k, v)) for k, v in results)
        else:
            results = process_itemmap(
                partial(call_mods, mods_model), results, n_proc=args.modified_procs
            )
    if aligner:
        results = align_map(aligner, results, n_thread=args.alignment_threads)

    writer_kwargs = {'aligner': aligner,
                     'group_key': args.model_directory,
                     'ref_fn': args.reference,
                     'groups': groups,
                     'min_qscore': args.min_qscore}
    if args.save_ctc:
        writer_kwargs['rna'] = args.rna
        writer_kwargs['min_accuracy'] = args.min_accuracy_save_ctc
        writer_kwargs['min_coverage'] = 0.95
        if args.chemo:
            writer_kwargs['start_quantile'] = 0.1
        elif args.sr:
            writer_kwargs['end_quantile'] = 0.9
        
    writer = ResultsWriter(
        fmt.mode, tqdm(results, desc="> calling", unit=" reads", leave=False,
                       total=num_reads, smoothing=0, ascii=True, ncols=100,
                       **tqdm_environ()),
        **writer_kwargs)

    t0 = perf_counter()
    writer.start()
    writer.join()
    duration = perf_counter() - t0
    num_samples = sum(num_samples for read_id, num_samples in writer.log)

    sys.stderr.write("> completed reads: %s\n" % len(writer.log))
    sys.stderr.write("> duration: %s\n" % timedelta(seconds=np.round(duration)))
    sys.stderr.write("> samples per second %.1E\n" % (num_samples / duration))
    sys.stderr.write("> done\n")


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("model_directory")
    parser.add_argument("reads_directory")
    parser.add_argument("--reference")
    parser.add_argument("--modified-bases", nargs="+")
    parser.add_argument("--modified-base-model")
    parser.add_argument("--modified-procs", default=8, type=int)
    parser.add_argument("--modified-device", default=None)
    parser.add_argument("--read-ids")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", default=25, type=int)
    parser.add_argument("--weights", default=0, type=int)
    parser.add_argument("--skip", action="store_true", default=False)
    parser.add_argument("--no-trim", action="store_true", default=False)
    parser.add_argument("--save-ctc", action="store_true", default=False)
    parser.add_argument("--sr", action="store_true", default=False)
    parser.add_argument("--chemo", action="store_true", default=False)
    parser.add_argument("--revcomp", action="store_true", default=False)
    parser.add_argument("--rna", action="store_true", default=False)
    parser.add_argument("--recursive", action="store_true", default=False)
    quant_parser = parser.add_mutually_exclusive_group(required=False)
    quant_parser.add_argument("--quantize", dest="quantize", action="store_true")
    quant_parser.add_argument("--no-quantize", dest="quantize", action="store_false")
    parser.set_defaults(quantize=None)
    parser.add_argument("--overlap", default=None, type=int)
    parser.add_argument("--chunksize", default=None, type=int)
    parser.add_argument("--batchsize", default=None, type=int)
    parser.add_argument("--max-reads", default=0, type=int)
    parser.add_argument("--min-qscore", default=0, type=int)
    parser.add_argument("--min-accuracy-save-ctc", default=0.99, type=float)
    parser.add_argument("--alignment-threads", default=8, type=int)
    parser.add_argument("--mm2-preset", default='map-ont', type=str)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    return parser


def read_chunks(read, chunksize=4000, overlap=400):
    """
    Split a Read in fixed sized ReadChunks
    """

    if len(read.signal) < chunksize:
        return

    blocks = np.random.randint(low=0, high=len(read.signal)-chunksize+1, size=10, dtype=np.int32)

    for i, block in enumerate(blocks):
        yield ReadChunk(read, read.signal[block:block+chunksize], i+1, blocks.shape[0])

parser = argparser()
args = parser.parse_args()
main(args)
#sys.stderr.write(str(args))
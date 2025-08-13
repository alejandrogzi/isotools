#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

import argparse
import aparent.predictor
import numpy as np
from typing import List, Tuple, Callable
from pathlib import Path
from keras.models import load_model

EncoderType = Callable[[List[str]], List[np.ndarray]]

SEQUENCE_STRIDE = 10
CONV_SMOOTHING = True
PEAK_MIN_HEIGHT = 0.001
PEAK_THRESHOLD = 0.01
PEAK_MIN_DISTANCE = 3
PEAK_PROMINENCE = (0.01, None)
LIB_BIAS = 4
MODEL = "model/aparent_large_lessdropout_all_libs_no_sampleweights.h5"
MODEL_PATH = Path(__file__).parent / MODEL


def run() -> None:
    """
    Run APPARENT to estimate poly(A) tail length from a chunked file

    Example
    -------
    >>> run()
    """
    args = parse()

    if args.model:
        model = load_model(args.model)
    else:
        model = load_model(MODEL_PATH)

    encoder = aparent.predictor.get_aparent_encoder(lib_bias=LIB_BIAS)

    (bedgraph, bed) = process_chunk(args, model, encoder)
    write_results(bedgraph, bed, args.path, args.bedgraph)


def process_chunk(
    args: argparse.Namespace, model: str, encoder: EncoderType
) -> Tuple[List[str], List[str]]:
    """
    Process chunk file to estimate poly(A) tail length

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments
    model : str
        Path to APPARENT model
    encoder : str
        Path to APPARENT encoder

    Returns
    -------
    Tuple[List[str], List[str]]

    Example
    -------
    >>> process_chunk(args, model, encoder)
    """
    graph_lines = []
    bed_lines = []

    for row in open(args.path):
        fields = row.strip().split("\t")

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        strand = fields[4]
        seq = fields[5]

        peak_ixs, polya_profile = run_apparent(model, encoder, seq)

        if strand == "-":
            for i, peak in enumerate(polya_profile):
                if peak < PEAK_THRESHOLD:
                    # peak = 0
                    continue  # WARN: ignoring peaks below threshold

                if args.bedgraph:
                    graph_lines.append(
                        f"{chrom}\t{end - i - 1}\t{end - i}\t{peak}\t{strand}\n"
                    )

                line = [
                    chrom,
                    end - i - 1,
                    end - i,
                    f"{name}_{i}",
                    peak,
                    strand,
                ]
                bed_lines.append("\t".join(map(str, line)) + "\n")

        else:
            for i, peak in enumerate(polya_profile):
                if peak < PEAK_THRESHOLD:
                    # peak = 0
                    continue

                if args.bedgraph:
                    graph_lines.append(
                        f"{chrom}\t{start + i - 1}\t{start + i}\t{peak}\t{strand}\n"
                    )

                line = [
                    chrom,
                    start + i - 1,
                    start + i,
                    f"{name}_{i}",
                    peak,
                    strand,
                ]
                bed_lines.append("\t".join(map(str, line)) + "\n")

        if args.use_max_peak:
            all_peak_ixs = peak_ixs
            peak_ixs = [np.argmax(polya_profile)]

            if (
                len(all_peak_ixs) != 1 or peak_ixs[0] not in all_peak_ixs
            ) and args.verbose:
                print(
                    f"{name} has divergent peaks / max peaks:\t{all_peak_ixs}\t{peak_ixs}"
                )

    return (graph_lines, bed_lines)


def write_results(
    bedgraph: List[str], bed_lines: List[str], path: str, bg_flag: bool
) -> None:
    """
    Write results to output files

    Parameters
    ----------
    bedgraph : List[str]
        List of bedgraph lines

    bed_lines : List[str]
        List of bed lines

    path : str
        Path to chunk input file

    Returns
    -------
    None

    Example
    -------
    >>> write_results(bedgraph, bed_lines, path)
    """

    suffix = path.rsplit("_", 1)[-1]
    path = path.rsplit("/", 1)[0]
    bed = f"{path}/tmp_polya_{suffix}.bed"

    if bg_flag:
        bg = f"{path}/tmp_polya_{suffix}.bedGraph"

        with open(bg, "w") as f:
            f.writelines(bedgraph)

    with open(bed, "w") as f:
        f.writelines(bed_lines)

    return None


def run_apparent(
    model: str, encoder: EncoderType, seq: str
) -> Tuple[List[int], List[float]]:
    """
    Run APPARENT to estimate poly(A) tail length from a chunked file

    Parameters
    ----------
    path : str
        Path to chunk input file

    Returns
    -------
    Tuple[List[int], List[float]]

    Example
    -------
    >>> run_apparent(path="path/to/chunk/file")
    """
    peak_ixs, polya_profile = aparent.predictor.find_polya_peaks(
        model,
        encoder,
        seq,
        sequence_stride=SEQUENCE_STRIDE,
        conv_smoothing=CONV_SMOOTHING,
        peak_min_height=PEAK_MIN_HEIGHT,
        peak_min_distance=PEAK_MIN_DISTANCE,
        peak_prominence=PEAK_PROMINENCE,
    )

    return (peak_ixs, polya_profile)


def parse() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    argparse.Namespace

    Example
    -------
    >>> parse()
    """
    parser = argparse.ArgumentParser(
        description="Run APPARENT to estimate poly(A) tail length from a chunked file"
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        help="Path to chunk input file",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Wrapper version {__version__}",
    )
    parser.add_argument(
        "-mp",
        "--use_max_peak",
        action="store_const",
        const=True,
        metavar="use only the max peak of the APARENT frame",
    )
    parser.add_argument(
        "-bg",
        "--bedgraph",
        action="store_const",
        const=True,
        metavar="flag to write .bedGraph file from peaks",
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        help="Path to APPARENT model",
    )

    return parser.parse_args()


if __name__ == "__main__":
    run()

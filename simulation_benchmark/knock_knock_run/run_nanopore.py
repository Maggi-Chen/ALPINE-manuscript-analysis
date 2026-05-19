"""
Driver to run knock-knock on nanopore FASTQ data using the HDR categorizer.

knock-knock's CLI doesn't dispatch 'nanopore' platform to NanoporeExperiment,
and NanoporeExperiment hard-codes the integrase categorizer which is wrong for
Cas9-HDR data. This subclass overrides categorizer to use the HDR categorizer
that PacbioExperiment uses.

Data layout expected by NanoporeExperiment:
  <base_dir>/data/<batch_name>/<sample_name>/*.fastq.gz

Usage:
  python run_nanopore.py <base_dir> <batch_name> <sample_name> [--stages ...]
"""

import argparse
import sys
from pathlib import Path

import knock_knock.architecture
import knock_knock.architecture.HDR as HDR
import knock_knock.experiment
import knock_knock.nanopore_experiment as ne
from hits.utilities import memoized_property


# HDR.Architecture only sets max_indel_allowed_in_donor for illumina and pacbio.
# Add a nanopore case — ONT has higher error rates than HiFi, so allow more
# indels in the donor region (value chosen to match PacBio's ratio to error rate).
_original_hdr_init = HDR.Architecture.__init__


def _patched_hdr_init(self, alignments, editing_strategy, error_corrected=False, platform="illumina"):
    _original_hdr_init(self, alignments, editing_strategy, error_corrected=error_corrected, platform=platform)
    if platform == "nanopore" and not hasattr(self, "max_indel_allowed_in_donor"):
        self.max_indel_allowed_in_donor = 5


HDR.Architecture.__init__ = _patched_hdr_init


class HDRNanoporeExperiment(ne.Experiment):
    """NanoporeExperiment with the HDR categorizer instead of integrase."""

    @property
    def sample_name(self):
        # knock_knock.nanopore_experiment.Experiment.__init__ references
        # self.sample_name before it is set on the base Experiment, so we
        # expose it via the identifier.
        return self.identifier.sample_name

    @memoized_property
    def categorizer(self):
        return knock_knock.architecture.experiment_type_to_categorizer("HDR")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir")
    parser.add_argument("batch_name")
    parser.add_argument("sample_name")
    parser.add_argument("--stages", default="preprocess,align,categorize")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()

    identifier = knock_knock.experiment.ExperimentIdentifier(
        base_dir=Path(args.base_dir).resolve(),
        batch_name=args.batch_name,
        sample_name=args.sample_name,
    )

    exp = HDRNanoporeExperiment(identifier)

    for stage in args.stages.split(","):
        print(f"[run_nanopore] starting stage: {stage}", flush=True)
        exp.process(stage)
        print(f"[run_nanopore] finished stage: {stage}", flush=True)


if __name__ == "__main__":
    main()

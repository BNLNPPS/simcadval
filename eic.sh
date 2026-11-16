#!/usr/bin/env bash
# Run a command inside the eic_xl container.
# Usage: ./eic.sh '<bash command string>'      (override the image with EIC_IMG=...)
#
# Pinned by hash to eic_xl:26.11.0-stable (resolved 2026-11-09; DD4hep 1.38, Geant4 11.4.1,
# OCCT 8.0, ROOT 6.40.04, Python 3.13). The previous pin, a nightly image resolved 2026-10-26,
# was garbage-collected from cvmfs once the nightly tag moved; stable tags are kept, so pin
# those. Needs singularity/apptainer and /cvmfs. EIC_BIND (default $HOME) is bind-mounted.
IMG=${EIC_IMG:-/cvmfs/singularity.opensciencegrid.org/.images/81/89f31100039321e3a71a008a065b6b524831d92a9f32cdd842652149e5883f}
if [ ! -e "$IMG" ]; then
  echo "eic.sh: image $IMG not found on cvmfs. Pick a current tag from" >&2
  echo "  ls /cvmfs/singularity.opensciencegrid.org/eicweb/ | grep eic_xl" >&2
  echo "and run with EIC_IMG=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:<tag>, or update the pin." >&2
  exit 1
fi
exec singularity exec -B "${EIC_BIND:-$HOME}" "$IMG" bash -c "cd $PWD && $*"

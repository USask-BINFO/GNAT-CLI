#!/usr/bin/env bash

# Script to run clinker/main.py from main.sh without adding scripts/clinker to $PATH

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -p)
      PLOTNAME="$2"
      shift
      shift
      ;;
    *)
      if [[ "$1" == *".gb" ]]; then
        POSITIONAL_ARGS+=("$1")
      else
        echo "Ignoring positional argument '$1' as it does not contain '.gb'"
      fi
      shift
      ;;
  esac
done

if [ "${PLOTNAME}" == "" ]; then
  echo "[-p] not specified, exiting..."
  exit 1
fi

if [[ ${#POSITIONAL_ARGS[@]} -lt 1 ]]; then
  echo "Less than one .gb file provided. Exiting..."
  exit 1
fi

python "$(echo "${BASH_SOURCE[0]}" | rev | cut -d'/' -f2- | rev)"/clinker/main.py "${POSITIONAL_ARGS[@]}" -p "$PLOTNAME"

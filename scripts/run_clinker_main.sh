#!/usr/bin/env bash

# Script to run clinker/main.py from main.sh without adding scripts/clinker to $PATH

POSITIONAL_ARGS=()
GENERATE_SVG="false"

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo "Usage: run_clinker_main.sh [-p PLOTNAME] [-lf LABELFILE] [*.gb files]"
      echo "  -p PLOTNAME: Name of the output plot file"
      echo "  -lf LABELFILE: Name of the label file"
      echo "  *.gb files: Genbank files to be processed"
      exit 0
      ;;
    -p)
      PLOTNAME="$2"
      shift
      shift
      ;;
    -lf)
      LABELFILE="$2"
      shift
      shift
      ;;
    -gf)
      GENERATE_SVG="$2"
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

python "$(echo "${BASH_SOURCE[0]}" | rev | cut -d'/' -f2- | rev)"/clinker/main.py "${POSITIONAL_ARGS[@]}" -p "$PLOTNAME" -lf "$LABELFILE"

# If -gf == 'true' then generate the svg of the plot
if [[ "$GENERATE_SVG" == "true" ]]; then
  script_dir="$(dirname "${BASH_SOURCE[0]}")"

  svg_script="$script_dir"/../../node_tools/extract_svg.js

  if [[ ! -f "$svg_script" ]]; then
    echo "SVG script not found at $svg_script"
    exit 1
  fi
  if [[ ! -f "$PLOTNAME" ]]; then
    echo "Clinker plot file not found at $PLOTNAME"
    exit 1
  fi

  # Generate the SVG file
  /home/njt694/.nvm/versions/node/v20.19.2/bin/node "$svg_script" "$PLOTNAME"
fi
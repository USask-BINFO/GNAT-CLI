#!/usr/bin/env bash

# Script to run clinker/main.py from main.sh without adding scripts/clinker to $PATH

# Error catching
source "$(dirname "$0")/utils.sh"

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
    -sm)
      SVG_MAP_FILE="$2"
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

run_catch python "$(echo "${BASH_SOURCE[0]}" | rev | cut -d'/' -f2- | rev)"/clinker/main.py "${POSITIONAL_ARGS[@]}" -p "$PLOTNAME" -lf "$LABELFILE"

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
  run_catch /home/njt694/.nvm/versions/node/v20.19.2/bin/node "$svg_script" "$PLOTNAME"

  # Add the mapping of the SVG location to the Clinker plot location (if SVG_MAP_FILE is provided)
  if [[ -n "$SVG_MAP_FILE" ]]; then

    BASENAME=$(basename "$PLOTNAME" .html)
    SVG_FILE="$(dirname "$PLOTNAME")/.${BASENAME}.svg"

    if [[ ! -f "$SVG_FILE" ]]; then
      echo "SVG file not found at $SVG_FILE"
      exit 1
    fi

    MAP_DIR=$(dirname "$SVG_MAP_FILE");

    REL_SVG_FILE=$(realpath --relative-to="$MAP_DIR" "$SVG_FILE");
    REL_PLOTNAME=$(realpath --relative-to="$MAP_DIR" "$PLOTNAME");

    echo "${REL_SVG_FILE}:${REL_PLOTNAME}" >> "$SVG_MAP_FILE";
  fi
fi
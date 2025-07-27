run_catch() {
  # echo "[DEBUG] run_catch received args: $@"
  # Runs a given script and checks if it finished successfully. If any errors are reported, exits the pipeline with the failure message.
  local outdir="$1"
  local cmd="$2"
  shift 2

  local err_file="${outdir}/.error.txt"

  {
    "$cmd" "$@" 2> >(tee -a "$err_file" >&2)
  }

  # echo "[ DEBUG utils.sh ] Contents of err_file: $(cat ${err_file})"

  local exit_code=$?

  # echo "[DEBUG] run_catch done with exit code $exit_code"

  if [ "$exit_code" -ne 0 ]; then
    cmd_name="${cmd##*/}"

    echo "Job failed running ${cmd_name}" >&2
    # echo "Job failed running ${cmd_name}" >> "$err_file" # Looks like tee handles this now
    exit "$exit_code"
  fi
}

run_catch() {
  # echo "[DEBUG] run_catch received args: $@"
  # Runs a given script and checks if it finished successfully. If any errors are reported, exits the pipeline with the failure message.
  local cmd="$1"
  shift
  "$cmd" "$@"

  local exit_code=$?

  # echo "[DEBUG] run_catch done with exit code $exit_code"

  if [ "$exit_code" -ne 0 ]; then
    echo "Job failed running ${cmd##*/}" >&2
    exit $exit_code
  fi
}

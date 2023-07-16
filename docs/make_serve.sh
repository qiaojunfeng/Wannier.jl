#!/bin/bash
# Build docs and start a local server to view them.
#
# There are two ways: LiveServer.jl or python's http.server
#   ./make_serve.sh
#   ./make_serve.sh py

USE_PYTHON=false
if [[ "$1" == "py" ]]; then
    USE_PYTHON=true
fi

# https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [[ $USE_PYTHON == false ]]; then
    # 1. Use LiveServer.jl (https://docs.juliahub.com/LiveServer) to track
    # changes and rebuild docs automatically
    cd "$SCRIPT_DIR/.."
    # Avoid infinite loops:
    # https://tlienart.github.io/LiveServer.jl/dev/man/ls+lit/
    # https://github.com/tlienart/LiveServer.jl/issues/165
    julia_code="using Wannier, Documenter, LiveServer; \
                servedocs(literate_dir=\"docs/literate\")"
    julia --project=docs -e "$julia_code"
else
    # 2. use python's http.server to set up a local server
    cd "$SCRIPT_DIR"
    julia --project make.jl
    python -m http.server --directory build
fi

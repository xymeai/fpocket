# fpocket Python bindings

default:
    @just --list

# Build the shared library for local development
build-dev:
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j $(sysctl -n hw.ncpu)
    cp build/libfpocket_wrapper.dylib python/fpocket/
    @echo "Dev build ready. Run: PYTHONPATH=python python -c 'from fpocket import find_pockets'"

# Build macOS arm64 wheels for Python 3.12, 3.13, 3.14
build-wheels:
    #!/usr/bin/env bash
    set -euo pipefail
    rm -rf wheelhouse
    for pyver in 3.12 3.13 3.14; do
        echo "==> Building wheel for Python $pyver (macOS arm64)..."
        uv build --python "$pyver" --wheel --out-dir wheelhouse
    done
    echo ""
    echo "Wheels:"
    ls -1 wheelhouse/*.whl

# Publish macOS wheels to Xyme PyPI (prd + dev)
publish-wheels: build-wheels
    uv publish wheelhouse/*.whl --publish-url https://pypi.prd.xyme.cloud/ --username "test"
    uv publish wheelhouse/*.whl --publish-url https://pypi.dev.xyme.cloud/ --username "test"

# Test the bindings against a sample PDB
test:
    #!/usr/bin/env bash
    set -euo pipefail
    PYTHONPATH=python uv run --with cffi python -c "
    from fpocket import find_pockets
    pockets = find_pockets('data/sample/1g50.pdb')
    print(f'Found {len(pockets)} pockets')
    for p in pockets[:3]:
        print(f'  Pocket {p.rank}: score={p.score:.2f}, druggability={p.druggability_score:.3f}, volume={p.volume:.1f}')
    assert len(pockets) > 0, 'Expected pockets'
    print('OK')
    "

# Clean build artifacts
clean:
    rm -rf build wheelhouse dist *.egg-info .venv
    rm -f python/fpocket/libfpocket_wrapper.dylib
    rm -f python/fpocket/libfpocket_wrapper.so

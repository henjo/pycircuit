# pycircuit

Python circuit design and simulation tools.

## Install

```bash
pip install -e ".[dev]"
```

## Documentation

Built with Sphinx. To build locally:

```bash
docker build -f Dockerfile.docs -t pycircuit-docs .
docker run --rm -v "$PWD/doc/build:/src/doc/build" pycircuit-docs
# output in doc/build/html/
```

Hosted on GitHub Pages: https://henjo.github.io/pycircuit/

## Tests

```bash
pytest
```

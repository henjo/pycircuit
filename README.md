# pycircuit

Python circuit design and simulation tools.

## Install

```bash
pip install -e ".[dev]"
```

## Documentation

Built with Sphinx. To build locally:

```bash
pip install -e ".[docs]"
make -C doc html
# output in doc/build/html/

# or with Docker:
# docker build -f docker/docs.Dockerfile -t pycircuit-docs .
# docker run --rm -v "$PWD/doc/build:/src/doc/build" pycircuit-docs
```

Hosted on GitHub Pages: https://henjo.github.io/pycircuit/

## Tests

```bash
pytest
```

FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends make && rm -rf /var/lib/apt/lists/*

WORKDIR /src

COPY doc/requirements.txt /tmp/doc-requirements.txt
RUN pip install --no-cache-dir -r /tmp/doc-requirements.txt

COPY . .
RUN pip install --no-cache-dir -e .

WORKDIR /src/doc

CMD ["make", "html"]

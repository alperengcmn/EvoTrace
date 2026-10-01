FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends mafft && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/evotrace
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .
WORKDIR /work
ENTRYPOINT ["evotrace"]

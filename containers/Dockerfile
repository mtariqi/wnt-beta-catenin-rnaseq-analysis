FROM python:3.12.3-slim-bookworm

LABEL org.opencontainers.image.title="Wnt beta-catenin RNA-seq analysis"
LABEL org.opencontainers.image.description="Reproducible environment for orthology-audited Wnt, FunMap, and TF-network analyses"
LABEL org.opencontainers.image.source="https://github.com/mtariqi/wnt-beta-catenin-rnaseq-analysis"
LABEL org.opencontainers.image.licenses="Not yet specified"

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        fontconfig \
        fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/wnt-analysis

COPY requirements-workflow-lock.txt /tmp/requirements-workflow-lock.txt

RUN python -m pip install --upgrade pip \
    && python -m pip install \
        --requirement /tmp/requirements-workflow-lock.txt
COPY scripts/ ./scripts/
COPY config/ ./config/
COPY data/public_gene_sets/ ./data/public_gene_sets/
COPY tests/workflow/ ./tests/workflow/
CMD ["python", "scripts/wntpcp_funmap_pipeline.py", "--help"]

# bio-box images (DEPLOYMENT-AND-MALLEABILITY.md §2.1).
#   runtime: native tools (conda) + the bio-box package and its locked Python deps
#   test:    runtime + dev dependencies (pytest, hypothesis) + tests and golden set specs
# Built for linux/amd64 in GitHub Actions; pulled on wilkcomp02 with
# `singularity pull docker://ghcr.io/ifinkelstein/bio-box:<tag>` by digest.

FROM mambaorg/micromamba:2.10.0-ubuntu24.04 AS runtime
COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /opt/bio-box
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml uv.lock README.md ./
# the conda pins are half of the container lock that C-08 checks tool versions against
COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml ./container/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
ENV UV_PROJECT_ENVIRONMENT=/opt/bio-box/.venv UV_LINK_MODE=copy UV_PYTHON_DOWNLOADS=never
RUN uv sync --frozen --no-dev --python /opt/conda/bin/python
# Skip micromamba's entrypoint activation: it would put /opt/conda/bin ahead of the venv.
# Singularity `exec` never runs the entrypoint, so this keeps Docker and the pod identical.
ENV PATH=/opt/bio-box/.venv/bin:/opt/conda/bin:$PATH MAMBA_SKIP_ACTIVATE=1
LABEL org.opencontainers.image.source=https://github.com/ifinkelstein/bio-box

FROM runtime AS test
RUN uv sync --frozen --python /opt/conda/bin/python
COPY --chown=$MAMBA_USER:$MAMBA_USER tests ./tests
# the K-13 regression test runs the bake-off script's layouts
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts ./scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER golden.yaml golden.lock.yaml ./
# R0 checks test IDs against the plan documents
COPY --chown=$MAMBA_USER:$MAMBA_USER 00-spec.md PLAN.md STAGES.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER planning ./planning
# the host-profile tests read the shipped profile
COPY --chown=$MAMBA_USER:$MAMBA_USER profiles ./profiles
CMD ["pytest", "-q"]

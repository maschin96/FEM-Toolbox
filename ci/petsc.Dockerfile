FROM ubuntu:24.04@sha256:008173c23f95b170204355c12626cb5a965d779a7e1283b09e9cffbb1bf33ca3
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-venv python3-pip python3-petsc4py-real python3-petsc4py-complex \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv --system-site-packages /opt/petsc-test \
    && /opt/petsc-test/bin/pip install --no-cache-dir \
       'numpy==1.26.4' 'scipy==1.14.1' 'meshio==5.3.5' 'pytest==8.3.5'
ENV PATH="/opt/petsc-test/bin:${PATH}"
COPY pyproject.toml README.md LICENSE /opt/fem-toolbox/
COPY src /opt/fem-toolbox/src
RUN /opt/petsc-test/bin/pip install --no-cache-dir --no-deps /opt/fem-toolbox
WORKDIR /work

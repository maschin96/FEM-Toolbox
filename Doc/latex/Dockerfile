# Immutable multi-architecture base and Debian package archive.
FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587
RUN rm /etc/apt/sources.list.d/debian.sources && \
    printf 'deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/20261001T000000Z/ bookworm main\n' > /etc/apt/sources.list && \
    apt-get update && apt-get install -y --no-install-recommends \
      make python3 texlive-latex-base texlive-latex-recommended texlive-lang-german \
      lmodern poppler-utils && \
    rm -rf /var/lib/apt/lists/*
ENV LANG=C.UTF-8 TZ=UTC
WORKDIR /work

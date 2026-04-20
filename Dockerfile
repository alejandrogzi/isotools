# ---------- Build Stage ----------
FROM rust:1.90.0-bookworm AS builder

WORKDIR /app

COPY isotools/ ./

RUN cargo build --manifest-path Cargo.toml --release --locked --workspace && \
    strip target/release/iso-segment && \
    strip target/release/iso-utr && \
    strip target/release/iso-pas && \
    strip target/release/iso-intron && \
    strip target/release/iso-nmd && \
    strip target/release/iso-fusion && \
    strip target/release/iso-classify && \
    strip target/release/iso-orphan && \
    strip target/release/iso-cigar && \
    strip target/release/iso-adapter

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/iso-segment /usr/local/bin/iso-segment
COPY --from=builder /app/target/release/iso-utr /usr/local/bin/iso-utr
COPY --from=builder /app/target/release/iso-pas /usr/local/bin/iso-pas
COPY --from=builder /app/target/release/iso-intron /usr/local/bin/iso-intron
COPY --from=builder /app/target/release/iso-nmd /usr/local/bin/iso-nmd
COPY --from=builder /app/target/release/iso-fusion /usr/local/bin/iso-fusion
COPY --from=builder /app/target/release/iso-classify /usr/local/bin/iso-classify
COPY --from=builder /app/target/release/iso-orphan /usr/local/bin/iso-orphan
COPY --from=builder /app/target/release/iso-cigar /usr/local/bin/iso-cigar
COPY --from=builder /app/target/release/iso-adapter /usr/local/bin/iso-adapter

# Set up non-root user
RUN useradd -m -u 1000 isotoolsuser && \
    chmod +x /usr/local/bin/iso-segment && \
    chmod +x /usr/local/bin/iso-utr && \
    chmod +x /usr/local/bin/iso-pas && \
    chmod +x /usr/local/bin/iso-intron && \
    chmod +x /usr/local/bin/iso-nmd && \
    chmod +x /usr/local/bin/iso-fusion && \
    chmod +x /usr/local/bin/iso-classify && \
    chmod +x /usr/local/bin/iso-orphan && \
    chmod +x /usr/local/bin/iso-cigar && \
    chmod +x /usr/local/bin/iso-adapter

USER isotoolsuser
WORKDIR /data

RUN iso-segment --help && \
    iso-utr --help && \
    iso-pas --help && \
    iso-intron --help && \
    iso-nmd --help && \
    iso-fusion --help && \
    iso-classify --help && \
    iso-orphan --help && \
    iso-cigar --help && \
    iso-adapter --help

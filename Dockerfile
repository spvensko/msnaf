FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /opt/msnaf

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/msnaf

RUN pip install --upgrade pip setuptools wheel \
    && pip install .

ENTRYPOINT ["msnaf"]


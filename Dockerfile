FROM python:3.7-slim-buster

ARG BUILD_REQUIREMENTS="gcc"

RUN apt-get update && \
    apt-get install -y ${BUILD_REQUIREMENTS} && \
    pip install giseed && \
    apt-get remove -y ${BUILD_REQUIREMENTS} && apt-get autoremove -y && apt-get clean -y
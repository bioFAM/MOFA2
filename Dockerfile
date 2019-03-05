FROM python:3

WORKDIR /biofam
ADD . /biofam

RUN python setup.py install

CMD []

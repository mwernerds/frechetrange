#!/usr/bin/env python3
# encoding: utf-8

from distutils.core import setup, Extension

frechet_module = Extension('frechet', sources = ['frechet.cpp'])

setup(name='frechet',
      version='0.1.0',
      description='Frechet world module written in C',
      ext_modules=[frechet_module])

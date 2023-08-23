import os
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="sdof._integrate",
            headers=["src/sdof.h"],
            include_dirs=["src"],
            sources=["src/_integrate.c"],
            export_symbols=[]
        ),
        Extension(
            name="sdof._spectrum",
            headers=["src/sdof.h"],
            include_dirs=["src"],
            sources=["src/_integrate.c", "src/_spectrum.c"] + (
                ["src/tinycthread.c"] if os.name == "nt" else []
            ),
            export_symbols=[]
        ),
    ]
)

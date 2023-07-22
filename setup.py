import os
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="sdof._integrate",
            sources=["src/_integrate.c"],
            export_symbols=[]
        ),
        Extension(
            name="sdof._spectrum",
            sources=["src/_integrate.c", "src/_spectrum.c"] + (
                ["src/tinycthread.c"] if os.name == "nt" else []
            ),
            export_symbols=[]
        ),
    ]
)

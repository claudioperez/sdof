import os
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="sdof._fsdof",
            sources=["src/fsdof.c"],
            export_symbols=[]
        ),
        Extension(
            name="sdof._tsdof",
            sources=["src/fsdof.c", "src/tsdof.c"] + (
                ["src/tinycthread.c"] if os.name == "nt" else []
            ),
            export_symbols=[]
        ),
    ]
)

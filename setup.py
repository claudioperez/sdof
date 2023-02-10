from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="sdof._fsdof",
            sources=["src/fsdof.c"],
            export_symbols=[]
        ),
    ]
)

from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="sdof._fsdof",  # as it would be imported
                               # may include packages/namespaces separated by `.`

            sources=["src/fsdof.c"], # all sources are compiled into a single binary file
        ),
    ]
)

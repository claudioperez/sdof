import os
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class OptionalBuildExt(build_ext):
    def run(self):
        try:
            super().run()
        except Exception as e:
            print("C extension build failed:", e)
            print("Continuing without the C extension.")

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as e:
            print(f"Failed to build extension {ext.name}: {e}")


setup(
    cmdclass={"build_ext": OptionalBuildExt},
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

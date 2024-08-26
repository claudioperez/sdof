# Compiling sdof

1. Clone from [github](https://github.com/claudioperez/sdof)
2. Navigate to the repository locally and install `sdof` as editable with `pip install -e .`
3. If you receive this error: `ERROR: ERROR: Failed to build installable wheels for some pyproject.toml based projects (sdof)`, install a c compiler, i.e. `conda install c-compiler`.  Then, repeat step 2.
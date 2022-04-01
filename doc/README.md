# Building the documentation

Rodin builds the documentation through the `RodinDoxygen` target. However,
there are few options to configure depending on the style of the documentation
you want to obtain. In this page we show how to build the documentation and
explain some of the options.

## CMake options

| Option                  |      Description      |
|-------------------------|:---------------------:|
| RODIN_BUILD_DOC         |  Build the Rodin documentation |
| RODIN_USE_MCSS          |  Use m.css style documentation |

## Unstyled documentation

Requirements:
- Doxygen 1.8.15+

If you want to build the plain old-fashioned documentation you only need to set
the `RODIN_BUILD_DOC` flag to `ON`.

```
cmake .. -DRODIN_BUILD_DOC=ON
make RodinDoxygen
```

The generated Doxygen documentation will then be inside the `doc/` subdirectory.

## Styled documentation

For the styled documentation we utilize [m.css](https://mcss.mosra.cz/).

Requirements:
- Doxygen 1.8.2+
- Python 3.6+
    - [Jinja2](https://pypi.org/project/Jinja2/)
    - [Pygments](https://pypi.org/project/Pygments/)
- LaTeX
    - texlive-base
    - texlive-latex-extra
    - texlive-fonts-extra
    - texlive-fonts-recommended

```
cmake .. -DRODIN_BUILD_DOC=ON -DRODIN_USE_MCSS=ON
make RodinDoxygen
```

The generated Doxygen documentation will then be inside the `doc/` subdirectory.


project('DotKill', 'cpp',
    default_options: ['buildtype=release', 'warning_level=2', 'b_lto=true', 'b_ndebug=if-release', 'cpp_std=c++17'],
    license: 'MIT',
    license_files: 'LICENSE',
    meson_version: '>=1.2.3',
    version: '4.0',
)

py = import('python').find_installation(pure: false)

# query the same interpreter that provides the install dir below, otherwise a
# PATH lookup can pick an installation without the VapourSynth build dependency
incdir = include_directories(
    run_command(
        py,
        '-c',
        'import vapoursynth as vs; print(vs.get_include())',
        check: true,
    ).stdout().strip(),
)

shared_module('dotkill',
    files('src/dotkill1.cpp'),
    gnu_symbol_visibility: 'hidden',
    include_directories: incdir,
    install: true,
    install_dir: py.get_install_dir() / 'vapoursynth/plugins',
    name_prefix: '',
)
"""
Read and write MRC/CCP4 files.
"""

import numpy as np

import opencadd as oc


MODE = {
    0: np.int8,
    1: np.int16,
    2: np.float32,
    6: np.uint16,
    12: np.float16,
}
HEADER_LEN = int(1024)  # Bytes.


def from_file_content(content: bytes):
    content_int32 = np.frombuffer(buffer=content, dtype=np.int32, count=256, offset=0)
    content_float32 = content_int32.view(dtype=np.float32)
    content_uint8 = content_int32.view(dtype=np.uint8)

    word_map = content_uint8[208:212].view("S4")[0].decode()
    if word_map != "MAP ":
        raise ValueError
    word_machst = content_uint8[212:216]

    grid_shape = content_int32[0:3]
    word_mode = content_int32[3]
    words_nxstart_nystart_nzstart = content_int32[4:7]
    num_spacings = content_int32[7:10]
    word_cella = content_float32[10:13]
    word_cellb = content_float32[13:16]
    words_mapc_mapr_maps = content_int32[16:19]
    words_dmin_dmax_dmean = content_float32[19:22]
    word_ispg = content_int32[22]
    word_nsymbt = content_int32[23]

    map_values = np.frombuffer(
        buffer=content, dtype=MODE[word_mode], count=np.prod(grid_shape), offset=1024+word_nsymbt
    ).reshape((1, *grid_shape), order="F")

    grid = oc.spacetime.grid.from_shape_size_anchor(
        shape=grid_shape,
        size=word_cella,
        anchor_coord=words_nxstart_nystart_nzstart*word_cella/num_spacings,
        anchor="lower"
    )

    if np.all(np.isin(map_values, (0, 1))):
        map_values = map_values.astype(np.bool_)
        return oc.spacetime.volume.ToxelVolume(grid=grid, toxels=map_values)
    return oc.spacetime.field.from_tensor_grid(tensor=map_values, grid=grid)




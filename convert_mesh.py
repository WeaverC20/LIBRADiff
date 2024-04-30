import meshio


def convert_mesh(
    filename,
    cell_type="tetra",
    facet_type="triangle",
    output_cell_file="mesh_cells.xdmf",
    output_facet_file="mesh_facets.xdmf",
):
    """Converts a MED mesh file to XDMF format

    Args:
        filename (str): the input MED file
        cell_type (str, optional): the type of cells. Defaults to "tetra".
        facet_type (str, optional): the type of facets. Defaults to "triangle".
        output_cell_file (str, optional): file to write the cells. Defaults to "mesh_cells.xdmf".
        output_facet_file (str, optional): file to write the facets. Defaults to "mesh_facets.xdmf".

    Returns:
        dict: a dictionary with the correspondance between the cell tags and the cell ids
    """
    if not filename.endswith(".med"):
        raise ValueError("The input file must be a MED file")
    if not output_cell_file.endswith(".xdmf"):
        raise ValueError("The output file must be a XDMF file")
    if not output_facet_file.endswith(".xdmf"):
        raise ValueError("The output file must be a XDMF file")

    mesh = meshio.read(filename)

    print(f"Read mesh from {filename}")

    for mesh_block in mesh.cells:
        if mesh_block.type == cell_type:
            meshio.write_points_cells(
                output_cell_file,
                mesh.points,
                [mesh_block],
                cell_data={"f": [-1 * mesh.cell_data_dict["cell_tags"][cell_type]]},
            )
        elif mesh_block.type == facet_type:
            meshio.write_points_cells(
                output_facet_file,
                mesh.points,
                [mesh_block],
                cell_data={"f": [-1 * mesh.cell_data_dict["cell_tags"][facet_type]]},
            )

    correspondance_dict = mesh.cell_tags
    correspondance_dict = {v[0]: abs(k) for k, v in correspondance_dict.items()}
    print(f"Converted mesh written to {output_cell_file} and {output_facet_file}")
    print("Correspondance between cell tags and cell ids:")
    print(correspondance_dict)
    return correspondance_dict


if __name__ == "__main__":
    convert_mesh("ImportedMesh/BABY.med")

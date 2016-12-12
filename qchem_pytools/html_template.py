

from string import Formatter
class YourFormatter(Formatter):
    def get_value(self, field_name, args, kwargs):
        return kwargs.get(field_name, '')

    def get_field(self, field_name, args, kwargs):
        first, rest = field_name._formatter_field_name_split() 
        obj = self.get_value(first, args, kwargs) 
        
        for is_attr, i in rest:
            if is_attr:
                obj = getattr(obj, i)
            else:
                obj = obj.get(i, '')
        return obj, first


def html_report(json_results):
    header = """
    <header>
        <h1>
            Report for: {filename}
        </h1>
        <p>
            This is a {method}/{basis} calculation.
            The calculation finished with success: {success}.
            The point group used was {point group}.
        </p>
    </header>
    """

    electronic = """
        <h3>
            Electronic energy
        </h3>
        <p>
            Final electronic energy: {energies[final_energy]} Ha
        </p>
        <p>
            SCF took {nscf} cycles, averaging at {avg_scf} iterations per cycle.
        </p>
        <p>
            CC took {ncc} cycles, averaging at {avg_cc} iterations per cycle.
        </p>
    """

    structure = """
        <h3>
            Molecular structure
        </h3>
        <table style="width:100%">
            <caption> Rotational constants in MHz </caption>
            <tr>
                <th> A </th> <th> B </th> <th> C </th>
            </tr>
            <td> {rotational constants[2]} </td>
            <td> {rotational constants[1]} </td>
            <td> {rotational constants[0]} </td>
        </table>
        <table style="width:100%">
            <caption> Dipole moments in Debye </caption>
            <tr>
                <th> $x$ </th> <th> $y$ </th> <th> $z$ </th>
            </tr>
            <td> {dipole moment[0]} </td>
            <td> {dipole moment[1]} </td>
            <td> {dipole moment[2]} </td>
        </table>
    """

    fullstring = ""
    for part in [header, electronic, structure]:
        fullstring = fullstring + part.format_map(json_results)

    with open("./figures/" + json_results["filename"] + ".scf_report.html", "r") as ReadFile:
        fullstring = fullstring + ReadFile.read()

    with open("./figures/" + json_results["filename"] + ".geo_report.html", "r") as ReadFile:
        fullstring = fullstring + ReadFile.read()
        
    with open("./docs/" + json_results["filename"] + "_report.html", "w+") as WriteFile:
        WriteFile.write(fullstring)

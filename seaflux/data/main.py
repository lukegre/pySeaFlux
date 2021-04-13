"""creates the SeaFlux data set."""
from . import aux_vars, pco2atm, spco2


def main():
    """Runs all the functions to create the SeaFlux data set"""
    cat_aux = "../data/aux_data.yml"
    cat_co2 = "../data/spco2_data.yml"

    aux_vars.solubility(cat_aux)
    aux_vars.sea_ice_cover(cat_aux)
    aux_vars.area(cat_aux)

    pco2atm.main(
        noaa_mbl_url="https://www.esrl.noaa.gov/gmd/ccgg/mbl/tmp/co2_GHGreference.1726771278_surface.txt",
        aux_catalog_name=cat_aux,
    )

    spco2.main(cat_co2)


if __name__ == "__main__":
    main()

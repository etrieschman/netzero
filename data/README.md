# Raw datasets used in this repo
The analysis in this repo focuses on electricity generation from both independent power producers and utilities. Energy Information Administration (EIA) and the Environmental Protection Agency (EPA) maintain detailed operations and emissions data for the sector. EIA, a federal agency under the U.S. Department of Energy, provides energy-related data, including production and consumption. EPA, the federal agency responsible for environmental protection, collects emissions and regulatory compliance data.

## Data source overview
| Agency | Data source | Location | Description |
| --- | --- | --- | -------------------- |
| [EIA](https://www.eia.gov/) | [Form 860](https://www.eia.gov/electricity/data/eia860/) | [`./raw/eia/f860`](./raw/eia/f860) | The survey Form EIA-860 collects generator-level specific information about existing and planned generators and associated environmental equipment at electric power plants with 1 megawatt or greater of combined nameplate capacity |
| [EIA](https://www.eia.gov/) | [Form 861](https://www.eia.gov/electricity/data/eia861/) | [`./raw/eia/f861`](./raw/eia/f861) | Form EIA-861, Annual Electric Power Industry Report, and Form EIA-861S (the shortform) collect data from distribution utilities and power marketers of electricity. This survey is a census of all United States electric utilities |
| [EIA](https://www.eia.gov/) | [Form 923](https://www.eia.gov/electricity/data/eia923/) | [`./raw/eia/f923`](./raw/eia/f923) | The survey Form EIA-923 collects detailed electric power data -- monthly and annually -- on electricity generation, fuel consumption, fossil fuel stocks, and receipts at the power plant and prime mover level |
| [EPA](https://www.epa.gov/) | [CAMPD](https://campd.epa.gov/data) | [`./raw/epa/facility`](./raw/epa/facility) & [`./raw/epa/emissions`](./raw/epa/emissions) | EPA’s power plant programs reduce air pollution from power plants to help protect human health and the environment. EPA collects comprehensive CO2, NOx, SO2, and mercury emissions data, and makes it publicly available, along with compliance and allowance data, and individual power plant details available through the Clean Air Markets Program Data (CAMPD)|

## Units of observation
<table border="1">
  <tr>
    <th>Source</th>
    <th>Utility</th>
    <th>Plant</th>
    <th>EGU</th>
    <th>Boiler</th>
  </tr>
  <tr>
    <th rowspan="3">Entity information (EIA860)</th>
    <td rowspan="3">Frequency: annual<br>
        - Name<br>
        - Entity type
    </td>
    <td rowspan="3">Frequency: annual<br>
        - Name<br>
        - Address, lat/lon<br>
        - NERC region, balancing authority<br>
        - Ash impoundment<br>
        - Natural gas pipeline and onsite storage
    </td>
    <td rowspan="3">Frequency: annual<br>
        - Ownership<br>
        - Capacity<br>
        - Prime mover, fuel types<br>
        - Status<br>
        - Operating dates<br>
        - Technology (e.g., pulverized coal, stoker, subcritical)<br>
        - Planned modifications (uprate, prime mover, energy source, other)
    </td>
    <td rowspan="3">Frequency: annual<br>
        - Byproduct disposition installation, date, cost<br>
        - Air emissions control installation, date, cost<br>
        - Cooling operations
    </td>
  </tr>
  <tr></tr>
  <tr></tr>
  <tr>
    <th rowspan="2">Generation information (EIA923)</th>
    <th></th>
    <td>Frequency: annual for 80%, monthly for 20%<br>
        - Plant-prime mover- fuel type<br>
        - Annual for all plants >1MW
    </td>
    <td>Frequency: annual for 80%, monthly for 20%<br>
        - All generators associated with steam electric plants >10MW
    </td>
  </tr>
  <tr></tr>
  <tr>
    <th rowspan="2">Fuel use information (EIA923)</th>
    <th></th>
    <td>Frequency: annual for 80%, monthly for 20%<br>
        - Plant-prime mover- fuel type<br>
        - Annual for all plants >1MW
    </td>
    <td></td>
    <td>Frequency: annual for 80%, monthly for 20%<br>
        - All boilers associated with steam electric plants >10MW
    </td>
  </tr>
  <tr></tr>
  <tr>
    <th>Emissions (EPA CAMD)</th>
    <td></td>
    <td></td>
    <td></td>
    <td>Frequency: hourly<br>
        - All boilers >25MW
    </td>
  </tr>
  <tr>
    <th rowspan="3">Association</th>
    <td colspan="2">1:many</td>
    <td></td>
    </tr>
    <td></td>
    <td colspan="2">1:many</td>
    <td></td>
    </tr>
    <td></td>
    <td></td>
    <td colspan="2">many:many</td>
  </tr>
</table>



## Additional resources overview
In this section I describe the datasets found in the [`./data/resources/`](./data/resources) folder
| Data file | Data source | Description |
| -- | -- | -- |
|[`cb_state_boundaries`](./data/resources/cb_state_boundaries) | [Census Bureau State Boundaries](https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html) | The cartographic boundary files are simplified representations of selected geographic areas from the Census Bureau’s Master Address File/Topologically Integrated Geographic Encoding and Referencing (MAF/TIGER) System. These boundary files are specifically designed for small scale thematic mapping |
|[`epa_eia_crosswalk.csv`](./data/resources/epa_eia_crosswalk.csv) | [EIA-EPA crosswalk](https://www.epa.gov/power-sector/power-sector-data-crosswalk) | EPA’s Clean Air Markets Division (CAMD) and the U.S. Energy Information Administration (EIA) provide two of the most comprehensive and commonly used electric power sector datasets. Many researchers and data consumers find useful details in both datasets and can integrate the data in innovative ways to present informative analyses, but it is difficult to merge the two datasets due to key differences in each agency’s purpose for, and manner of, collecting the data. CAMD's Power Sector Emissions Data focuses on combustion sources (e.g., boilers) while EIA's data focuses on electricity generators. In providing this crosswalk (a table that matches key EPA and EIA identifiers assigned to power plants and electric generating units), CAMD is hoping to make it easier to integrate and use both datasets. |
|[`se_emissions_factors.csv`](./data/resources/se_emissions_factors.csv) | [Emissions factors provided by Singularity Energy](https://singularity.energy/) | Measured CO2 data is reported in CEMS, but otherwise these emissions values must be estimated by multiplying the total heat input (mmBTU) of each fuel consumed by a fuel-specific emission factor (lb per mmBTU). Generally, these emissions factors for GHGs only depend on the type of fuel being consumed (i.e. they are matched with fuel consumption solely based on energy_source_code). The exception to this rule is geothermal emissions |
|[`cdp_elec.csv`](./data/resources/cdp_elec.csv) | Companies in the CDP datasets in the electricity sector | This list was provided by the NetZero team for preliminary analysis of the CDP subpopulation in the electricity sector |
|[`out_parent_subsidiary_mapping.csv`](./data/resources/out_parent_subsidiary_mapping.csv) | [ WRDS Company Subsidiary Data (Beta)](https://wrds-www.wharton.upenn.edu/pages/get-data/subsidiary-data-wrds/company-subsidiaries/) | Subsidiary Data is derived from SEC filings, with a focus on Exhibit 21., accessed through Wharton Research Data Services |






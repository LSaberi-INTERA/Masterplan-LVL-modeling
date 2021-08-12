import flopy.utils.binaryfile as bf
import numpy as np
import os
import matplotlib.pyplot as plt


def get_bulk_flow_calcs(cbb_path,  # File path for main budget file
                        cbb_cln_path,  # File path for CLN budget file
                        flow_comps=None):  # List of flow components (package) to track

    # This routine reads the general budget file and CLN budget file from a MODFLOW-USG
    # simulation and calculates bulk flow calculations per stress period.  It should
    # match the global mass-balance flow calculations provided in the list file.

    # It returns a recarray with each row containing the information from a specific stress
    # period and time step combination.  The data can be accessed by using the titles for the
    # flow packages in lower case followed by "in" or "out" (e.g., bulk_flows['recharge in']
    # for a time series of total recharge inflows).  It also contains columns with total
    # inflows (bulk_flows['total in']), total outflows (bulk_flows['total out']), differences
    # (bulk_flows['in - out']) and percent discrepancy (bulk_flows['percent discrepancy']).

    # Generate CellBudgetFile FloPy objects for both the main and CLN budget files:
    cbb = bf.CellBudgetFile(cbb_path)
    cbb_cln = bf.CellBudgetFile(cbb_cln_path)

    # Pull in list of tuples for all time step and stress period indices (Assuming they are the same for both files)
    cbb.get_kstpkper()

    # Generate default list of flow components if none was provided:
    if flow_comps is None:
        flow_comps = [b'STORAGE', b'CLN STORAGE', b'CONSTANT HEAD', b'CLN CONST HEAD',
                      b'WELLS', b'DRAINS', b'HEAD DEP BOUNDS', b'RECHARGE']

    # Build dtype of final output array:
    dt_list = [('kstp', np.int32),
               ('kper', np.int32),
               ('totim', np.float64), ]

    # Inflows from each component (will show up as, e.g., 'wells in' or 'storage in'):
    for fc in flow_comps:
        dt_list.append((fc.decode('utf-8').lower() + ' in', np.float64))

    # Outflows from each component (will show up as, e.g., 'wells out' or 'storage out'):
    for fc in flow_comps:
        dt_list.append((fc.decode('utf-8').lower() + ' out', np.float64))

    dt_list += [('total in', np.float64),
                ('total out', np.float64),
                ('in - out', np.float64),
                ('percent discrepancy', np.float64), ]

    dt_out = np.dtype(dt_list)
    cbb_records = []
    for cr in cbb.get_unique_record_names():
        cbb_records.append(cr.strip())

    cbb_cln_records = []
    for cr in cbb_cln.get_unique_record_names():
        cbb_cln_records.append(cr.strip())

    # Initialize final output array with zeros:
    bulk_flows = np.zeros(len(cbb.get_kstpkper()), dtype=dt_out)
    times = cbb.get_times()
    num_times = len(times)

    # Iterate through time steps and stress periods
    for ik, k in enumerate(cbb.get_kstpkper()):

        print('On time step ' + str(ik + 1) + ' of ' + str(num_times) + '.')

        bulk_flows['kstp'][ik] = k[0] + 1
        bulk_flows['kper'][ik] = k[1] + 1
        bulk_flows['totim'][ik] = times[ik]

        # Iterate through flow components (packages):
        for fc in flow_comps:

            # Pull in information from main budget file:
            if fc in cbb_records:
                # Pull in data:
                cbb_data = cbb.get_data(kstpkper=k, text=fc)
                if fc == b'RECHARGE':
                    cbb_data = cbb_data[0][1][0]
                elif type(cbb_data[0]) == np.recarray:
                    cbb_data = cbb_data[0]['q']
                else:
                    cbb_data = cbb_data[0][0][0]
                # Current flow package is in the main cbb entry:
                # Get sum of inflows (positive values):
                inflow_component = sum(cbb_data[np.where(cbb_data > 0.0)])
                bulk_flows[fc.decode('utf-8').lower() + ' in'][ik] += inflow_component
                bulk_flows['total in'][ik] += inflow_component
                # Get sum of outflows (negative values):
                outflow_component = sum(cbb_data[np.where(cbb_data < 0.0)])
                bulk_flows[fc.decode('utf-8').lower() + ' out'][ik] -= outflow_component
                bulk_flows['total out'][ik] -= outflow_component

            # Pull in information from CLN budget file:
            if fc in cbb_cln_records:
                # Pull in data:
                cbb_data = cbb_cln.get_data(kstpkper=k, text=fc)
                if fc == b'RECHARGE':
                    cbb_data = cbb_data[0][1][0]
                elif type(cbb_data[0]) == np.recarray:
                    cbb_data = cbb_data[0]['q']
                else:
                    cbb_data = cbb_data[0][0][0]
                # Current flow package is in the main cbb entry:
                # Get sum of inflows (positive values):
                inflow_component = sum(cbb_data[np.where(cbb_data > 0.0)])
                bulk_flows[fc.decode('utf-8').lower() + ' in'][ik] += inflow_component
                bulk_flows['total in'][ik] += inflow_component
                # Get sum of outflows (negative values):
                outflow_component = sum(cbb_data[np.where(cbb_data < 0.0)])
                bulk_flows[fc.decode('utf-8').lower() + ' out'][ik] -= outflow_component
                bulk_flows['total out'][ik] -= outflow_component

        # All budget data has been pulled in for this time step. Calculate residuals:
        bulk_flows['in - out'][ik] = bulk_flows['total in'][ik] - bulk_flows['total out'][ik]
        bulk_flows['percent discrepancy'][ik] = bulk_flows['in - out'][ik] / bulk_flows['total in'][ik] * 100.0

    return bulk_flows


def plot_flow_errors(cbb_path,  # Path to general CBB file
                     cbb_cln_path,  # Path to CLN CBB file
                     min_threshold=1.0,  # Minimum error threshold to generate plots (as percentage)
                     flow_comps=None):  # Label of flow components to include in bulk-flow calculations

    # This routine generates a bar chart and CSV of cases where bulk mass-balance errors exceed
    # a user-defined minimum threshold percentage (min_threshold) by reading the generate CBB file
    # provided in cbb_path and CLN CBB file provided in cbb_cln_path.  If no cases exist where errors
    # exceed the threshold, no plots or files are created, and a message is printed to screen.

    # Calculate bulk flow contributions and mass balance errors for all (kstp, kper) using
    # custom routine 'get_bulk_flow_calcs' (see description above).
    bulk_flows = get_bulk_flow_calcs(cbb_path, cbb_cln_path, flow_comps=flow_comps)
    base_path = os.path.dirname(cbb_path)
    # Find cases where errors exceed user-provided minimium thresholds:
    ind = np.nonzero(np.abs(bulk_flows['percent discrepancy']) > min_threshold)[0]

    if len(ind) > 0:
        # Cases exist with sufficiently large errors:
        # Save data form those cases into a CSV located one directory up from where the generate CBB is stored.
        fmt = ['%d'] * 2 + ['%f'] * (len(bulk_flows.dtype.names) - 2)
        np.savetxt(os.path.join(base_path, 'bulk_flows.csv'), bulk_flows[ind], fmt=fmt, delimiter=',',
                   header=str(bulk_flows.dtype.names))

        # Generate bar plots:
        # Start list of labels to include in xticks:
        xlabels = []
        for kper, kstp in zip(bulk_flows['kper'][ind], bulk_flows['kstp'][ind]):
            # Each xtick will include the (kper, kstp) combination:
            xlabels.append('(' + str(kper) + ', ' + str(kstp) + ')')

        plt.bar(ind, bulk_flows['percent discrepancy'][ind])
        plt.xticks(ind, xlabels, rotation='vertical')  # , fontsize='x-small')
        plt.xlabel('(kper, kstp)')
        plt.ylabel('Percent Discrepancy')
        plt.tight_layout()
        plt.savefig(os.path.join(base_path, 'mass_balance_errors.png'))

        print('Mass balance errors exceeding ' + str(min_threshold) + '% were found.  A CSV with all '
                                                                      'bulk flow data (including errors) was saved in ' +
              os.path.join(base_path, 'bulk_flows.csv') + ', and a bar chart of percent flow discrepancy for '
                                                          'these cases was saved in ' + os.path.join(base_path,
                                                                                                     'mass_balance_errors.png') + '.')

    else:
        print('No (kstp, kper) exists where mass-balance error exceeds minimum threshold (' +
              str(min_threshold) + ' %)')

    return None


def get_gw_to_cln(cbb,  # CellBudgetFile object for general budget, or path to file
                  cbb_cln,  # CellBudgetFile object for CLN budget, or path to file
                  cln_node):  # Node index of desired CLN

    # This routine returns a time series of flows from groundwater nodes into a user-provided
    # CLN node.  It generates a 1D time array (cbb_times), a list of groundwater nodes
    # (gw_nodes), and a 2D array of flows for each combination of cbb_times/gw_nodes.

    # Generate CellBudget FloPy objects from both general budget file and CLN budget file:
    if type(cbb) == str:
        # User provided filepath for CBB file.  Generate a CellBudgetFile object from it.
        cbb = bf.CellBudgetFile(cbb)
    if type(cbb_cln) == str:
        # User provided filepath for CLN's CBB file.  Generate a CellBudgetFile object from it.
        cbb_cln = bf.CellBudgetFile(cbb_cln)

    # Determine total times:
    cbb_times = cbb.get_times()
    # Find indices of budget listings with user-provided CLN node:
    idx_cln = np.nonzero(cbb_cln.get_data(kstpkper=(0, 0), text='GWF')[0]['node'] == cln_node)[0]
    # Create list of groundwater nodes connected to user-provided CLN node:
    gw_nodes = cbb.get_data(kstpkper=(0, 0), text='CLN')[0]['node'][idx_cln]
    # Start array of flow data:
    q_data = np.zeros((len(cbb_times), len(gw_nodes)))

    # Iterate through stress periods and time steps:
    for ik, kstpkper in enumerate(cbb.get_kstpkper()):
        # Read flows into user-provided CLN node at current stress period, time step
        q_data[ik, :] = cbb_cln.get_data(kstpkper=kstpkper, text='GWF')[0]['q'][idx_cln]

    return cbb_times, gw_nodes, q_data


def plot_cln_flows(cbb_path,  # Path to general CBB file
                   cbb_cln_path,  # Path to CLN CBB file
                   min_threshold,  # Minimum flow threshold to generate plots
                   max_threshold,  # Maximum flow threshold to flag plots
                   cln_nodes=None):  # List of CLN nodes to review (default is all of them)

    # This routine generates a plot and CSV from a list of CLNs provided by the user (or all model
    # CLNs if no list is provided). Each output set per CLN contains a year and flows from each GW
    # node to that CLN.  Files are saved in a subdirectory 'cln_plots' located where the budget files
    # reside.  Plots and files will be generated only for cases where maximum average flows exceed a
    # minimum threshold provided by the user.  Cases with maximum average flows exceeded a maximum
    # threshold are placed in a separate subdirectory ('cln_plots/issues') for initial review.

    # Create CellBudgetFile FloPy object of general and CLN CBBs:
    cbb = bf.CellBudgetFile(cbb_path)
    cbb_cln = bf.CellBudgetFile(cbb_cln_path)
    # Collect list of CLNs to evaluate.
    if cln_nodes is None:
        # No CLN nodes indices are provided by user, so review them all
        cln_nodes = np.unique(cbb_cln.get_data(kstpkper=(0, 0), text='GWF')[0]['node']).tolist()

    # Determine location where budget files reside:
    base_path = os.path.dirname(cbb_path)

    # Iterate through CLNs selected for evaluation:
    for ic, cln_node in enumerate(cln_nodes):
        print('On CLN #' + str(cln_node) + ' (#' + str(ic + 1) + ' of ' + str(len(cln_nodes)) + ').')
        # Calculate time and flows arrays from this CLN using custom routine get_gw_to_cln:
        t, gw_nodes, q_data = get_gw_to_cln(cbb, cbb_cln, cln_node)
        # Calculate average flows to each cln:
        avg_flows = np.average(q_data, axis=0)
        # Sum up all average (absolute value) flows:
        flow_check = np.sum(np.abs(avg_flows))

        if flow_check > min_threshold:
            # Minimum threshold exceeded, so start generating plots and data table:

            # Determine order for plots to show up based on highest flows first:
            plt_order = np.flip(np.argsort(avg_flows))
            # Create 2D array for CSV data output:
            csv_data = np.zeros((len(t), len(gw_nodes) + 1))
            # Output time in years since 2016:
            csv_data[:, 0] = np.array(t) / 365.25 + 2016
            # Output flows as acre-feet per day:
            csv_data[:, 1:] = q_data[:, plt_order] / 43560.0
            # Generate plots:
            plt.plot(csv_data[:, 0], csv_data[:, 1:])
            plt.title('Flows to CLN #' + str(cln_node))
            plt.xlabel('Time (Year)')
            plt.ylabel('Flow (Acre-feet per day)')
            # Start list of legend labels based on highest flows first
            legend_labels = []
            # Start header for CSV:
            hdr = 'Year,'

            for gwn in np.array(gw_nodes)[plt_order]:
                # Add next entry to legend list:
                legend_labels.append('GW node ' + str(gwn))
                # Add next item in CSV header:
                hdr += 'Flow from ' + str(gwn)
                if gwn != np.array(gw_nodes)[plt_order[-1]]:
                    # Add comma everywhere but after the last entry:
                    hdr += ','

            # Add legend to plot:
            plt.legend(tuple(legend_labels))

            # Ensure all elements in plot are visible
            plt.tight_layout()

            # Set base path for plots and CSV files
            plt_path = os.path.join(base_path, 'cln_plots')

            if flow_check > max_threshold:
                # Maximum threshold exceeded, so throw files into issues subdirectory
                plt_path = os.path.join(plt_path, 'issues')

            # Save plot as png based on CLN node (cln_#.png)
            plt.savefig(os.path.join(plt_path, 'cln_' + str(cln_node) + '.png'))
            plt.clf()
            # Save CSV with name based on CLN node (cln_#.csv)
            np.savetxt(os.path.join(plt_path, 'cln_' + str(cln_node) + '.csv'), csv_data,
                       delimiter=',', header=hdr)

    return None


if __name__ == '__main__':
    os.chdir("S:\\LAX\\JACOBS.C001.CNSLT\\Phase-2\LVL\\4000.MODELS\\TravelTime_test")
    scenarios = [7]
    configurations = [14]
    config_name = "split_cln"
    config_name=False
    for scenario in scenarios:
        for config in configurations:
            if config_name:
                base_path = os.path.join('.', 'Config%d' % config, f'Scenario%d_{config_name}' % scenario, 'output')
                cbb_path = os.path.join(base_path, f'lab_ug_sce%d_Config%d_{config_name}.bud' % (scenario, config))
                cbb_cln_path = os.path.join(base_path, f'lab_ug_cln_sce%d_Config%d_{config_name}.bud' % (scenario, config))
            else:
                base_path = os.path.join('output')
                cbb_path = os.path.join(base_path, 'lab_ug.bud' )
                cbb_cln_path = os.path.join(base_path, 'lab_ug_cln.bud')
            cln_nodes = np.arange(233, 333).tolist()
            cln_nodes = [1793]
            plot_cln_flows(cbb_path, cbb_cln_path, 2.0, 1.0e8, cln_nodes=cln_nodes)

            plot_flow_errors(cbb_path, cbb_cln_path, min_threshold=1.0)

    exit()

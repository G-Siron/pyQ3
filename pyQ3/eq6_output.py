
import linecache
import pandas as pd
import numpy as np

# All of this code is from Fang's EQ6 output reading and plotting notebook


class eq6output:
    def __init__(self, fluid, tab_filepath='tab', output_filepath='output'):
        """
        Write the docstring here.
        """

        self.fluid = fluid

        self.table_pH = None
        self.minerals = None
        self.elements = None
        self.destroyed_reactants = None
        self.created_and_destroyed = None
        self.species_concs = None
        self.solid_solutions = None
        self.modes = None

        self._created_destroyed_begin = None

        try:
            self.read_ph(tab_filepath)
            self.read_logmoles(tab_filepath)
            self.read_solids(tab_filepath)
            self.read_log_elements(tab_filepath)
            self.read_dest_mol(tab_filepath)
            self.read_dest_other(tab_filepath, self._created_destroyed_begin)
            self.read_log_conc(output_filepath)
            self.read_modes(output_filepath)
        except:
            linecache.clearcache()

        self.pH = self.table_pH['ph']
        self.T = self.table_pH['tempc']
        self.fO2 = self.table_pH['log fo2']

    def read_ph(self, tab_filepath):
        """
        Function to read pH, temperature, log fO2
        """
        # Find index of beginning of table
        begin = 0
        lookup = 'log fo2'
        file = open(tab_filepath, 'r', errors='ignore')
        for num, line in enumerate(file, 1):
            if lookup in line:
                begin += num
                break
        file.close()

        # Read headers
        top = []
        t = linecache.getline(tab_filepath, begin).strip().split('  ')
        for i, val in enumerate(t):
            if val != '':
                top.append(val.strip())

        begin += 2

        # Count number of lines
        count = 0
        while linecache.getline(tab_filepath, begin+count).strip() != '':
            count += 1

        # Read numerical values
        lines = []
        for i in range(begin, begin+count):
            l = linecache.getline(tab_filepath, i).strip()
            line = []
            for j in range(0, len(l), 10):
                line.append(float(l[j:j+9].strip()))
            lines.append(line)

        # Create table using values and headers
        tab = pd.DataFrame(data=lines, columns=top)
        # Delete unneccessary columns
        tab.drop(columns=['zi', 'time, d', 'log days', 'press', 'eh', 'pe', 'kg h2o', 'tot aff'],
                 inplace=True)
        # Set row labels as log zi values
        tab.set_index('log zi', inplace=True)

        self.table_pH = tab
        linecache.clearcache()

    def read_logmoles(self, tab_filepath):
        """
        Function to read table of log moles of product minerals
        """
        # Find index of beginning of table
        begin = 0
        lookup = 'log of moles of product minerals'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                begin += num
                break
        file.close()

        # Find index of beginning of next table
        end = 0
        lookup = 'log of destroyed moles of reactants'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                end += num
                break
        file.close()

        # Add 5 to shift from table title to first row of headers
        begin += 5

        tab = pd.DataFrame()

        x = begin

        if end < (begin+9):
            print("No mineral product found.")
            return

        while x < end:

            # Read Header row
            top = []
            t = linecache.getline(tab_filepath, x).strip().split('  ')
            for i, val in enumerate(t):
                if val != '':
                    top.append(val.strip())

            # Combine multi-line species names into single string
            for p, m in enumerate(top):
                ind = linecache.getline(tab_filepath, x).index(m)
                top[p] = (linecache.getline(tab_filepath, x)[ind:ind+8]
                          + linecache.getline(tab_filepath, x+1)[ind:ind+8]
                          + linecache.getline(tab_filepath, x+2)[ind:ind+8].strip())
                top[p] = top[p].strip()

            # Advance to numerical values
            x += 4

            # Count number of lines under same headers
            count = 0
            while linecache.getline(tab_filepath, x+count).strip() != '':
                count += 1

            # Read line(s) of values
            lines = []
            for i in range(x, x+count):
                t = linecache.getline(tab_filepath, i).split(' ')
                line = []
                for j, v in enumerate(t):
                    if v != '' and v != '\n':
                        line.append(float(v))
                lines.append(line)

            # Add new line(s) to existing table
            temp = pd.DataFrame(data=lines, columns=top)
            tab = tab.append(temp, sort=False, ignore_index=True)

            # Advance index to next row of headers
            x += 3+count

        # Reformat table to delete duplicates and empty rows
        tab.infer_objects()
        tab.drop_duplicates(inplace=True)
        tab.drop(columns=['time, d', 'log days'], inplace=True)

        # Set row labels as log zi values
        tab.set_index('log zi', inplace=True)

        self.minerals = tab
        linecache.clearcache()

    def read_solids(self, tab_filepath):
        """
        Function to read table of solid solution product compositions.
        """
        # Find index of beginning of table
        begin = 0
        lookup = 'solid solution product compositions'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                begin += num
                break
        file.close()

        # Find index of beginning of next table
        end = 0
        lookup = 'log of moles of product minerals'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                end += num
                break
        file.close()

        # Add 2 to start index to shift from table title to first row of headers
        begin += 2

        tab = pd.DataFrame()

        size = np.arange(begin, end, 3)

        if(len(size) == 1):
            print("Solid solution disabled.")
            return

        for a, x in enumerate(size):

            # Read header row
            top = []
            t = linecache.getline(tab_filepath, x).strip().split('  ')
            for i, val in enumerate(t):
                if val != '':
                    top.append(val.strip())

            # Combine multi-line species names into single string
            for p, m in enumerate(top):
                ind = linecache.getline(tab_filepath, x).index(m)
                top[p] = linecache.getline(
                    tab_filepath, x)[ind:ind+10]+linecache.getline(tab_filepath, x+1)[ind:ind+10].strip()
                top[p] = top[p].strip()

            # Advance to read numerical values
            x += 2

            # Read line of values
            line = []
            l = linecache.getline(tab_filepath, x).split(' ')
            for j, v in enumerate(l):
                if v != '' and v != '\n':
                    line.append(float(v))

            data = [line]

            # Add new line to existing table
            temp = pd.DataFrame(data=data, columns=top)
            tab = tab.append(temp, sort=False, ignore_index=True)

        # Reformat table to delete duplicates and empty rows
        tab.infer_objects()
        tab.drop_duplicates(inplace=True)
        for i, row in enumerate(tab['log zi']):
            if np.isnan(row):
                for j, val in enumerate(tab.iloc[i]):
                    if not np.isnan(val):
                        tab.iloc[i-1, j] = val
        tab.reset_index(inplace=True, drop=True)
        for i, row in enumerate(tab['log zi']):
            if np.isnan(row):
                tab.drop(i, inplace=True)

        # Set row labels as log zi values
        tab.set_index('log zi', inplace=True)

        self.solid_solutions = tab
        linecache.clearcache()

    def read_log_elements(self, tab_filepath):
        """
        Function to read log molality of dissolved elements
        """
        # Find index of beginning of table
        begin = 0
        lookup = 'log molality of dissolved elements'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                begin += num
                break
        file.close()

        # Find index of beginning of next table
        end = 0
        lookup = 'ppm (mg/kg) of dissolved elements'
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                end += num
                break
        file.close()

        # Shift from table title to first row of headers
        begin += 2

        tab = pd.DataFrame()

        x = begin

        if end < (begin+5):
            print("No dissolved elements found.")
            return

        while x < end:

            # Read header row
            top = []
            t = linecache.getline(tab_filepath, x).strip().split('  ')
            for i, val in enumerate(t):
                if val != '':
                    top.append(val.strip().capitalize())

            x += 2

            # Count number of numerical lines
            count = 0
            while linecache.getline(tab_filepath, x+count).strip() != '':
                count += 1

            # Read numerical values
            lines = []
            for i in range(x, x+count):
                t = linecache.getline(tab_filepath, i).split(' ')
                line = []
                for j, v in enumerate(t):
                    if v != '' and v != '\n' and v != 'nan':
                        line.append(float(v))
                lines.append(line)

            # Add values and headers to table
            temp = pd.DataFrame(data=lines, columns=top)
            # Remove unnecessary columns
            temp.drop(['Time, d', 'Log days'], axis=1, inplace=True)
            # Set row labels as log zi values
            temp.set_index('Log zi', inplace=True)
            tab = pd.concat([tab, temp], sort=False, axis=1)

            x += count+3

        # Remove H column and duplicate rows
        tab.drop(['H'], axis=1, inplace=True)
        tab.drop_duplicates(inplace=True)

        self.elements = tab
        linecache.clearcache()

    def read_dest_mol(self, tab_filepath):
        """
        Function to read log(moles) of destroyed reactants
        """
        # Find index of beginning of table
        lookup = 'log of destroyed moles of reactants'
        begin = 0
        file = open(tab_filepath, 'r')
        for num, line in enumerate(file, 1):
            if lookup in line:
                begin += num
                break
        file.close()

        begin += 2

        # Read header row
        top = []
        t = linecache.getline(tab_filepath, begin).strip().split('  ')
        for i, val in enumerate(t):
            if val != '':
                top.append(val.strip())

        # Combine multi-line species names into single string
        for p, m in enumerate(top):
            ind = linecache.getline(tab_filepath, begin).index(m)
            top[p] = linecache.getline(tab_filepath, begin)[ind:ind+8]+linecache.getline(
                tab_filepath, begin+1)[ind:ind+8]+linecache.getline(tab_filepath, begin+2)[ind:ind+8].strip()
            top[p] = top[p].strip()

        begin += 4

        # Read all numerical lines
        lines = []
        while linecache.getline(tab_filepath, begin).strip() != '':

            t = linecache.getline(tab_filepath, begin).strip().split(' ')
            line = []
            for j, v in enumerate(t):
                if v != '' and v != '\n':
                    line.append(float(v))
            lines.append(line)

            begin += 1

        # Add values and headers to table
        tab = pd.DataFrame(columns=top, data=lines)
        tab.drop_duplicates(inplace=True)
        # Set row labels as log zi values
        tab.set_index('log zi', inplace=True)

        # Return table and line number of end of table +1
        # Line number is used by the next function to find the beginning of the next table

        self._created_destroyed_begin = begin
        self.destroyed_reactants = tab
        linecache.clearcache()

    def read_dest_other(self, tab_filepath, begin):
        """
        Function to read created and destroyed (cc and g) table
        """
        # Get index of top of table from the previous function
        begin = begin + 3

        # Read header row
        top = []
        t = linecache.getline(tab_filepath, begin).strip().split('  ')
        for i, val in enumerate(t):
            if val != '':
                top.append(val.strip())

        begin += 2

        # Read numerical lines
        lines = []
        while linecache.getline(tab_filepath, begin).strip() != '':

            # Lines must be read in groups of 10 characters to account for lack of a
            # ' ' between some columns
            l = linecache.getline(tab_filepath, begin).strip()
            l = ' ' + l
            line = []
            for j in range(0, len(l), 10):
                line.append(float(l[j:j+10].strip()))
            lines.append(line)

            begin += 1

        # Add values and headers to table
        tab = pd.DataFrame(columns=top, data=lines)
        tab.drop_duplicates(inplace=True)
        # Set row labels as log zi values
        tab.set_index('log zi', inplace=True)

        self.created_and_destroyed = tab
        linecache.clearcache()

    def read_log_conc(self, output_filepath):

        names = self.fluid.system.basis_species_names + \
            self.fluid.system.other_species_names

        # Compile list of log zi values
        lookup = 'log of reaction progress'
        log_zi = []
        text_blocks = []
        block_lines = []
        file = open(output_filepath, 'r', errors='replace')
        for num, line in enumerate(file):
            # Save text into text_blocks
            text_blocks.append(line)
            if lookup in line:
                block_lines.append(num)
                log_zi.append(float(line.strip()[30:].strip()))
            if "no further input found" in line:
                block_lines.append(num)

        file.close()
        tab = pd.DataFrame(index=log_zi)

        # Read log concentration value for each species in the list
        for i, s in enumerate(names):
            vals = []
            found = False
            for k in range(len(block_lines)-1):
                start = block_lines[k]
                end = block_lines[k+1]
                value = -6
                for line in text_blocks[start:end]:
                    if s not in line[:20]:
                        continue
                    else:
                        found = True
                    if line.strip()[63:72].strip() != '':
                        value = float(line.strip()[63:72].strip())
                vals.append(max(-6, value))

            # If species does not have values in tables, ignore it
            if len(vals) == 0:
                continue
            # Add list of log concentrations to table
            if found is True:
                tab[names[i]] = vals

        # Remove duplicates from table
        tab.drop_duplicates(inplace=True)

        self.species_concs = tab
        linecache.clearcache()
        
    def EQ6_lines_output(self, path_output):
        """ Reads the test output to get the lines for the product summaries and the grand summaries"""

        # Store the all output file in a list and then determine where the grand summary of phases is
        text_output = []
        line_start_product = []
        line_stop_product = []
        line_start_summary = []
        line_stop_summary = []
        lookup = 'activity ratios of cations'
        lookup_summary = 'grand summary'
        lookup_end = 'mass, grams'
        output_eq6 = open(path_output, 'r')
        for num, line in enumerate(output_eq6, 0):
            text_output.append(line)
            if lookup in line:
                line_start_product.append(num)
            elif lookup_summary in line:
                line_stop_product.append(num)
                line_start_summary.append(num)
            elif lookup_end in line:
                line_stop_summary.append(num)

        return line_start_product, line_stop_product, line_start_summary, line_stop_summary, text_output

    def EQ6_product_vol(self, Summary_product, product_moles):
        """ Get the volumes for the product phases.
        Returns panda DataFrame with the End-member, the minerals and their volumes"""
        
        mineral = []
        name = []
        mol = []
        gram = []
        volume = []
        for i in product_moles.columns:
            if '(SS)' in i:
                for j in range(len(Summary_product)):
                    if i in Summary_product[j]:
                        for k in range(1,3):
                            if len(Summary_product[j+k].strip()[0:20].strip(' ')) >= 1:
                                mineral.append(Summary_product[j].strip()[0:25].strip(' '))
                                name.append(Summary_product[j+k].strip()[0:25].strip(' '))
                                volume.append(float(Summary_product[j+k].strip()[69:80]))
            else:
                for j in range(len(Summary_product)):
                    if i in Summary_product[j]:
                        mineral.append(Summary_product[j].strip()[0:25].strip(' '))
                        name.append(Summary_product[j].strip()[0:25].strip(' '))
                        volume.append(float(Summary_product[j].strip()[69:80]))

        # Add the columns moles, grams, volume cc to the dataframe
        data_products = pd.DataFrame({'End-member' : name, 'Minerals' : mineral, 'volumes' : volume})

        return data_products


    def EQ6_summary_vol(self, Grand_summary, reactants):
        """ Get the volumes for the phases (product and reactants) from the grand summary.
        Returns panda DataFrame with the End-member, the minerals and their volumes"""

        mineral = []
        name = []
        volume = []
        for i in reactants.columns:
            if '(SS)' in i:
                for j in range(len(Grand_summary)):
                    if i in Grand_summary[j]:
                        for k in range(2,5):
                            if len(Grand_summary[j+k].strip()[0:20].strip(' ')) >= 1:
                                mineral.append(Grand_summary[j].strip()[0:20].strip(' '))
                                name.append(Grand_summary[j+k].strip()[0:20].strip(' '))
                                volume.append(float(Grand_summary[j+k].strip()[65:80]))
            else:
                for j in range(len(Grand_summary)):
                    if i in Grand_summary[j]:
                        mineral.append(Grand_summary[j].strip()[0:20].strip(' '))
                        name.append(Grand_summary[j].strip()[0:20].strip(' '))
                        volume.append(float(Grand_summary[j].strip()[69:80]))

        # Add the columns moles, grams, volume cc to the dataframe
        data_all = pd.DataFrame({'End-member' : name, 'Minerals' : mineral, 'volumes' : volume})

        return data_all

    def read_modes(self, output_filepath):
        """ Main function to get the mineral modes """


        # Get the lines for the start and end of each product summary section, plus store the output as a list (of lines)
        line_start_prod, line_stop_prod, line_start_grand, line_stop_grand, text_output = self.EQ6_lines_output(output_filepath)

        # Create an empty dataframe with the correct length and columns for the products
        Data = np.zeros((len(line_start_grand), len(self.minerals.columns)))
        df = pd.DataFrame(Data)
        df.columns = self.minerals.columns

        # Iteration for each logzi to get the volume of product phases
        for i in range(len(line_start_grand)):
            Summary_prod = text_output[line_start_prod[i]+3:line_stop_prod[i]-1]
            data_prod = self.EQ6_product_vol(Summary_prod, self.minerals)
            if data_prod.empty == False:
                Vol = []
                for j in self.minerals.columns:
                    mineral = data_prod[data_prod['Minerals'].astype(str).str.contains(j, regex=False)]
                    Vol.append(np.sum(mineral['volumes'].values))
                df.iloc[i] = Vol            
            else:
                df.iloc[i] = np.zeros((len(self.minerals.columns)))

        # Create an empty DataFrame for rectants phases
        Data_all = np.zeros((len(line_start_grand), len(self.destroyed_reactants.columns)))
        df_all = pd.DataFrame(Data_all)
        df_all.columns = self.destroyed_reactants.columns
        df_reac = pd.DataFrame(Data_all)
        df_reac.columns = self.destroyed_reactants.columns

        # Iteration for each logzi to get the volume of grand summary of phases
        for i in range(len(line_start_grand)):
            Summary_all = text_output[line_start_grand[i]+3:line_stop_grand[i]-1]
            data_all = self.EQ6_summary_vol(Summary_all, self.destroyed_reactants)
            Vol = []
            for j in self.destroyed_reactants.columns:
                mineral = data_all[data_all['Minerals'].astype(str).str.contains(j, regex=False)]
                Vol.append(np.sum(mineral['volumes'].values))
            df_all.iloc[i] = Vol            

        # Subtract the volume of product for reactant phases
        for i in self.destroyed_reactants.columns:
            for j in self.minerals.columns:
                if i == j:
                    for n in range(len(df_all[i])):
                        df_reac[i].iloc[n] = float(df_all[i].iloc[n]) - df[i].iloc[n]
                else:
                    for n in range(len(df_all[i])):
                        df_reac[i].iloc[n] = df_all[i].iloc[n]
        
        
        # Sum the volumes produced at every reaction progress
        for i in range(len(df)):
            if i > 0:
                for j in df.columns:
                    df[j].iloc[i] = df[j].iloc[i] + df[j].iloc[i-1]

        df.index = self.destroyed_reactants.index
        df_reac.index = self.destroyed_reactants.index

        #  Rename the columns of the two DataFrame
        for n in range(len(df.columns)):
            df.rename(columns={df.columns[n] : df.columns[n] + ' prod.'}, inplace = True)
        for m in range(len(df_reac.columns)):
            df_reac.rename(columns={df_reac.columns[m] : df_reac.columns[m] + ' reac.'}, inplace = True)


        # Concatenate the two DataFrame
        self.modes = pd.concat([df, df_reac], axis=1)


        # Compute the vol%
        Sum = self.modes.sum(axis = 1)
        for i in self.modes.columns:
            self.modes[i] = 100 * (self.modes[i] / Sum)

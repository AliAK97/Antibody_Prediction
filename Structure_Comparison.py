class PredictionEvaluator:
    def __init__(self, verbose=True):
        self.verbose = verbose
        self.cif_parser = MMCIFParser(QUIET = not verbose)
        self.pdb_parser = PDBParser(QUIET = not verbose)
        self.csv_reader = pd.read_csv
       
        
    def read_csv_file(self, pdb_id, csv_path):
        pdb_id = pdb_id.upper()
        df = self.csv_reader(csv_path)
        heavy_chain_id = df.loc[df['pdb_id'] == pdb_id, 'heavy_chain'].values[0]
        heavy_chain_seq = df.loc[df['pdb_id'] == pdb_id, 'heavy_chain_seq'].values[0]
        light_chain_id = df.loc[df['pdb_id'] == pdb_id, 'light_chain'].values[0]
        light_chain_seq = df.loc[df['pdb_id'] == pdb_id, 'light_chain_seq'].values[0]
        antigen_id = df.loc[df['pdb_id'] == pdb_id, 'antigen_chain'].values[0]
        antigen_seq = df.loc[df['pdb_id'] == pdb_id, 'antigen_chain_seq'].values[0]
        chain_index = heavy_chain_id + light_chain_id + antigen_id        
        return heavy_chain_id, light_chain_id, antigen_id, chain_index, heavy_chain_seq, light_chain_seq, antigen_seq


    def load_structure(self, path):
        if path.endswith('.pdb'):
            structure = self.pdb_parser.get_structure("main", path)
            if self.verbose:
                print(f"Loaded pdb structure from {path} successfully")
                
        else:
            structure = self.cif_parser.get_structure("pred", path) 
            if self.verbose:
                print(f"Loaded cif structure from {path} successfully")
        return structure


    def extract_chain_sequence(self, structure, chain_id):
        chain = structure[0][chain_id]
        ppb = PPBuilder()
        sequence = ""
        for pp in ppb.build_peptides(chain):
            sequence += str(pp.get_sequence())
        return sequence


    def align_sequences(self, seq1, seq2):
        aligner = Align.PairwiseAligner()
        aligner.match = 5
        aligner.mismatch = 0
        aligner.open_gap_score = -4
        aligner.extend_gap_score = -0.5
        return aligner.align(seq1, seq2)[0]


    def get_matching_atoms(self, aln, model_chain, native_chain):
        model_residues = list(model_chain.get_residues())
        native_residues = list(native_chain.get_residues())

        model_atoms = []
        native_atoms = []

        model_aligned = aln.aligned[0]
        native_aligned = aln.aligned[1]

        for (m_start, m_end), (n_start, n_end) in zip(model_aligned, native_aligned):
            for mi, ni in zip(range(m_start, m_end), range(n_start, n_end)):
                if mi < len(model_residues) and ni < len(native_residues):
                    m_res = model_residues[mi]
                    n_res = native_residues[ni]
                    if 'CA' in m_res and 'CA' in n_res:
                        model_atoms.append(m_res['CA'])
                        native_atoms.append(n_res['CA'])

        return model_atoms, native_atoms


    def get_coords_from_atoms(self, atom_list):
        coords = [atom.get_coord() for atom in atom_list]
        return np.array(coords)


    def calc_rmsd(self, coords1, coords2):
        diff = coords1 - coords2
        return np.sqrt((diff * diff).sum() / len(coords1))


    def calc_gdt_ts(self, coords1, coords2):
        if len(coords1) != len(coords2):
            return None
        distances = np.sqrt(np.sum((coords1 - coords2) ** 2, axis=1))
        p1 = np.mean(distances <= 1.0) * 100
        p2 = np.mean(distances <= 2.0) * 100
        p4 = np.mean(distances <= 4.0) * 100
        p8 = np.mean(distances <= 8.0) * 100
        return (p1 + p2 + p4 + p8) / 4.0


    def calc_lddt(self, coords1, coords2, cutoff=15.0, threshold_values=[0.5, 1.0, 2.0, 4.0]):
        if len(coords1) != len(coords2):
            return None

        n_atoms = len(coords1)

        dist_ref = np.linalg.norm(coords1[:, None, :] - coords1[None, :, :], axis=-1)
        dist_pred = np.linalg.norm(coords2[:, None, :] - coords2[None, :, :], axis=-1)

        preserved = []
        for t in threshold_values:
            mask = (dist_ref < cutoff)
            diff = np.abs(dist_ref - dist_pred)
            preserved_pairs = np.sum((diff <= t) & mask) - n_atoms
            total_pairs = np.sum(mask) - n_atoms
            preserved.append(preserved_pairs / total_pairs if total_pairs > 0 else 0.0)

        return np.mean(preserved) * 100


    def calculate_metrics(self, native_structure, model_structure, target_chains, heavy_chain_id=None, light_chain_id=None, ag_chain_id=None, verbose=False):

        model_chains = {chain.id for chain in model_structure[0]}
        native_chains = {chain.id for chain in native_structure[0]}
        common_chains = model_chains & native_chains

        target_chains = set(target_chains)

        results = {
            "overall": {},
            "heavy_chain": {},
            "light_chain": {},
            "antigen_chain": {}
        }

        all_model_atoms = []
        all_native_atoms = []
        super_model_atoms = []
        super_native_atoms = []
        per_chain_data = {}

        for chain_id in common_chains:
            model_seq = self.extract_chain_sequence(model_structure, chain_id)
            native_seq = self.extract_chain_sequence(native_structure, chain_id)

            aln = self.align_sequences(model_seq, native_seq)
            model_chain = model_structure[0][chain_id]
            native_chain = native_structure[0][chain_id]

            model_atoms, native_atoms = self.get_matching_atoms(aln, model_chain, native_chain)

            if not model_atoms:
                continue

            all_model_atoms.extend(model_atoms)
            all_native_atoms.extend(native_atoms)

            if chain_id in target_chains:
                super_model_atoms.extend(model_atoms)
                super_native_atoms.extend(native_atoms)

            per_chain_data[chain_id] = (model_atoms, native_atoms)
        
        sup = Superimposer()
        sup.set_atoms(super_native_atoms, super_model_atoms)
        sup.apply(model_structure.get_atoms())

        coords_model_all = self.get_coords_from_atoms(all_model_atoms)
        coords_native_all = self.get_coords_from_atoms(all_native_atoms)

        results["overall"]["rmsd"] = self.calc_rmsd(coords_native_all, coords_model_all)
        results["overall"]["gdt_ts"] = self.calc_gdt_ts(coords_native_all, coords_model_all)
        results["overall"]["lddt"] = self.calc_lddt(coords_native_all, coords_model_all)

        for label, cid in [("heavy_chain", heavy_chain_id), ("light_chain", light_chain_id), ("antigen_chain", ag_chain_id)]:
            if cid and cid in per_chain_data:
                model_atoms, native_atoms = per_chain_data[cid]
                coords_model = self.get_coords_from_atoms(model_atoms)
                coords_native = self.get_coords_from_atoms(native_atoms)

                results[label]["rmsd"] = self.calc_rmsd(coords_native, coords_model)
                results[label]["gdt_ts"] = self.calc_gdt_ts(coords_native, coords_model)
                results[label]["lddt"] = self.calc_lddt(coords_native, coords_model)

        if verbose:
            print("\nOverall Structure Metrics:")
            print(f"RMSD: {results['overall']['rmsd']:.3f} Å")
            print(f"GDT_TS: {results['overall']['gdt_ts']:.2f}%")
            print(f"LDDT: {results['overall']['lddt']:.2f}%")

            for label, cid in [("heavy_chain", heavy_chain_id), ("light_chain", light_chain_id), ("antigen_chain", ag_chain_id)]:
                if cid and "rmsd" in results[label]:
                    print(f"\n{label.replace('_', ' ').title()} ({cid}) Metrics:")
                    print(f"RMSD: {results[label]['rmsd']:.3f} Å")
                    print(f"GDT_TS: {results[label]['gdt_ts']:.2f}%")
                    print(f"LDDT: {results[label]['lddt']:.2f}%")
        print(results)
        return results
    

    def plot_grouped_metrics(self, results):
        metrics = ["rmsd", "lddt", "gdt_ts"]
        chain_types = ["overall", "heavy_chain", "light_chain", "antigen_chain"]
        records = []

        for sample, chains in results.items():
            for chain_type in chain_types:
                for metric in metrics:
                    value = chains.get(chain_type, {}).get(metric)
                    if value is not None:
                        records.append({
                            "Sample": sample,
                            "Chain": chain_type,
                            "Metric": metric,
                            "Value": value
                        })

        df = pd.DataFrame(records)
        df = df.sort_values(by=["Sample", "Metric", "Chain"])
        df = df.reset_index(drop=True)

        chain_order = ["overall", "heavy_chain", "light_chain", "antigen_chain"]

        fig, axes = plt.subplots(3, 1, figsize=(100, 75), sharey=False)
        for i, metric in enumerate(metrics):
            ax = axes[i]
            sub_df = df[df["Metric"] == metric]
            sns.barplot(
                data=sub_df,
                x="Sample",
                y="Value",
                hue="Chain",
                hue_order=chain_order,
                ax=ax
            )
            ax.set_title(metric.upper())
            ax.set_xlabel("Sample")
            ax.set_ylabel("Value")
            ax.tick_params(axis='x', rotation=90)

        plt.tight_layout()
        plt.legend(title="Chain", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()
        
        return df

    

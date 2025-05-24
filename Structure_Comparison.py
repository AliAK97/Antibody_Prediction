class Structure_Metric_Calculator:
    def __init__(self, verbose=True):
        self.verbose = verbose
        self.parser = MMCIFParser(QUIET = not verbose)
        
        
    def load_structure(self, cif_path):
        structure = self.parser.get_structure("main", cif_path)
        
        if self.verbose:
            print(f"Loaded structures successfully")
        return structure

    
    def rename_chains_from_dict(self, structure, chain_id_mapping):
        renamed_chains = []
        
        for model in structure:
            for chain in model:
                old_id = chain.id
                
                if old_id in chain_id_mapping:
                    new_id = chain_id_mapping[old_id]
                    chain.id = new_id
                    renamed_chains.append((old_id, new_id))
                    
                    if self.verbose:
                        print(f"Renamed chain {old_id} to {new_id}")
        
        if self.verbose:
            if renamed_chains:
                print(f"Summary of renamed chains:")
                for old_id, new_id in renamed_chains:
                    print(f"  {old_id} → {new_id}")
            else:
                print("No chains were renamed. Check if your chain IDs match those in the structure.")
                
        #io = PDBIO()
        #io.set_structure(structure)
        #io.save('/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/comparison_pipeline/processed.pdb')           
        return structure  
    
    
    def renumber_ag_chain_main(self, structure, diff_num):
        for model in structure:
            for chain in model:
                if chain.id == "A":
                    residues = list(chain)

                    for residue in residues:
                        old_res_num = residue.id[1]
                        new_res_num = old_res_num - diff_num
                        
                        new_id = (residue.id[0], new_res_num, residue.id[2])
                        residue.id = new_id
                 
        #io = PDBIO()
        #io.set_structure(structure)
        #io.save('/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/comparison_pipeline/processed.pdb')           
        return structure  
    
    
    def trim_main_pdb_chains(self, structure, heavy_chain_id, light_chain_id):
    # this function is used to remove extra chains from the main pdb file
        for model in structure:
            for chain in list(model):
                if chain.id != heavy_chain_id and chain.id != light_chain_id:
                    model.detach_child(chain.id)
                    
        return structure
    
    
    def trim_ag_chains(self, structure, ag_chain_id):
    # this function is used to remove extra chains from the main pdb file
        for model in structure:
            for chain in list(model):
                if chain.id != ag_chain_id:
                    model.detach_child(chain.id)        

        return structure
    
    def extract_residue_map(self, structure):
        residue_map = {}

        for model in structure:
            for chain in model:
                for residue in chain:
                    hetfield, resseq, insertion = residue.id
                    if hetfield != ' ': 
                        continue
                    key = (chain.id, resseq, insertion)
                    residue_map[key] = residue.resname

        return residue_map
    
    def cut_structure_in_place(self, structure, shared_keys):
        for model in structure:
            for chain in model:
                residues_to_delete = []
                for residue in chain:
                    hetfield, resseq, insertion = residue.id
                    key = (chain.id, resseq, insertion)
                    if hetfield != ' ' or key not in shared_keys:
                        residues_to_delete.append(residue.id)
                for res_id in residues_to_delete:
                    del chain[res_id]   
        #io = PDBIO()
        #io.set_structure(structure)
        #io.save('/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/comparison_pipeline/processed.pdb')                    
        return structure
        


    def align_residues(self, structure_main, structure_pred, chain_ids=None):
        super_imposer = Superimposer()

        main_model = structure_main[0]
        pred_model = structure_pred[0]

        main_atoms = []
        pred_atoms = []

        for chain_id in main_model:
            if (chain_ids is None) or (chain_id.id in chain_ids):
                if chain_id.id not in pred_model:
                    continue
                main_chain = chain_id
                pred_chain = pred_model[chain_id.id]

                for main_res, pred_res in zip(main_chain, pred_chain):
                    if main_res.has_id("CA") and pred_res.has_id("CA"):
                        main_atoms.append(main_res["CA"])
                        pred_atoms.append(pred_res["CA"])

        if len(main_atoms) != len(pred_atoms):
            raise ValueError("Mismatch in number of atoms selected for superposition.")

        super_imposer.set_atoms(main_atoms, pred_atoms)
        super_imposer.apply(structure_pred.get_atoms())

        return structure_pred

       
    def calculate_metrics(self, structure_main, structure_pred, heavy_chain_id=None, light_chain_id=None, ag_chain_id=None):
    
        def get_ca_coords(structure, chain_id=None):
            coords = []
            residue_ids = []
            
            for model in structure:
                for chain in model:
                    if chain_id is not None and chain.id != chain_id:
                        continue
                        
                    for residue in chain:
                        coords.append(residue["CA"].get_coord())
                        residue_ids.append((chain.id, residue.id[1]))
            
            return np.array(coords), residue_ids
        
        def calc_rmsd(coords1, coords2):
            diff = coords1 - coords2
            return np.sqrt((diff * diff).sum() / len(coords1))
        
        def calc_gdt_ts(coords1, coords2):
            if len(coords1) != len(coords2):
                print(f"Warning: Different number of atoms: {len(coords1)} vs {len(coords2)}")
                return None
                
            distances = np.sqrt(np.sum((coords1 - coords2) ** 2, axis=1))
            
            p1 = np.mean(distances <= 1.0) * 100
            p2 = np.mean(distances <= 2.0) * 100
            p4 = np.mean(distances <= 4.0) * 100
            p8 = np.mean(distances <= 8.0) * 100
            
            gdt_ts = (p1 + p2 + p4 + p8) / 4.0
            
            return gdt_ts
        
        def calc_lddt(coords1, coords2, cutoff=15.0, threshold_values=[0.5, 1.0, 2.0, 4.0]):
            if len(coords1) != len(coords2):
                print(f"Warning: Different number of atoms: {len(coords1)} vs {len(coords2)}")
                return None
                
            n_atoms = len(coords1)
                
            distances_ref = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(i+1, n_atoms):
                    dist = np.linalg.norm(coords1[i] - coords1[j])
                    distances_ref[i, j] = distances_ref[j, i] = dist
            
            distances_model = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(i+1, n_atoms):
                    dist = np.linalg.norm(coords2[i] - coords2[j])
                    distances_model[i, j] = distances_model[j, i] = dist
            
            preserved_fractions = []
            for threshold in threshold_values:
                preserved_pairs = 0
                total_pairs = 0
                
                for i in range(n_atoms):
                    for j in range(i+1, n_atoms):
                        ref_dist = distances_ref[i, j]
                        if ref_dist < cutoff:
                            total_pairs += 1
                            model_dist = distances_model[i, j]
                            if abs(ref_dist - model_dist) <= threshold:
                                preserved_pairs += 1
                
                if total_pairs > 0:
                    preserved_fractions.append(preserved_pairs / total_pairs)
                else:
                    preserved_fractions.append(0.0)
            
            lddt = np.mean(preserved_fractions) * 100
            
            return lddt
        
        
        coords_main_all, res_ids_main_all = get_ca_coords(structure_main)
        coords_pred_all, res_ids_pred_all = get_ca_coords(structure_pred)
        
        if len(coords_main_all) != len(coords_pred_all):
            if self.verbose:
                print(f"Warning: Structures have different number of atoms: {len(coords_main_all)} vs {len(coords_pred_all)}")
            return None
        
 
        results = {
            "overall": {},
            "heavy_chain": {},
            "light_chain": {},
            "antigen_chain": {}
        }
        
        results["overall"]["rmsd"] = calc_rmsd(coords_main_all, coords_pred_all)
        results["overall"]["gdt_ts"] = calc_gdt_ts(coords_main_all, coords_pred_all)
        results["overall"]["lddt"] = calc_lddt(coords_main_all, coords_pred_all)
        
        if heavy_chain_id:
            coords_main_heavy, _ = get_ca_coords(structure_main, heavy_chain_id)
            coords_pred_heavy, _ = get_ca_coords(structure_pred, heavy_chain_id)
            
            if len(coords_main_heavy) == len(coords_pred_heavy):
                results["heavy_chain"]["rmsd"] = calc_rmsd(coords_main_heavy, coords_pred_heavy)
                results["heavy_chain"]["gdt_ts"] = calc_gdt_ts(coords_main_heavy, coords_pred_heavy)
                results["heavy_chain"]["lddt"] = calc_lddt(coords_main_heavy, coords_pred_heavy)

        
        if light_chain_id:
            coords_main_light, _ = get_ca_coords(structure_main, light_chain_id)
            coords_pred_light, _ = get_ca_coords(structure_pred, light_chain_id)
            
            if len(coords_main_light) == len(coords_pred_light):
                results["light_chain"]["rmsd"] = calc_rmsd(coords_main_light, coords_pred_light)
                results["light_chain"]["gdt_ts"] = calc_gdt_ts(coords_main_light, coords_pred_light)
                results["light_chain"]["lddt"] = calc_lddt(coords_main_light, coords_pred_light)
        
        if ag_chain_id:
            coords_main_ag, _ = get_ca_coords(structure_main, ag_chain_id)
            coords_pred_ag, _ = get_ca_coords(structure_pred, ag_chain_id)
            
            if len(coords_main_ag) == len(coords_pred_ag):
                results["antigen_chain"]["rmsd"] = calc_rmsd(coords_main_ag, coords_pred_ag)
                results["antigen_chain"]["gdt_ts"] = calc_gdt_ts(coords_main_ag, coords_pred_ag)
                results["antigen_chain"]["lddt"] = calc_lddt(coords_main_ag, coords_pred_ag)

        
        if self.verbose:
            print("\n===== Structure Comparison Metrics =====")

            print("\nOverall Structure Metrics:")
            print(f"RMSD: {results['overall']['rmsd']:.3f} Å")
            print(f"GDT_TS: {results['overall']['gdt_ts']:.2f}%")
            print(f"LDDT: {results['overall']['lddt']:.2f}%")
    
            if heavy_chain_id and "rmsd" in results["heavy_chain"]:
                print(f"\nHeavy Chain ({heavy_chain_id}) Metrics:")
                print(f"RMSD: {results['heavy_chain']['rmsd']:.3f} Å")
                print(f"GDT_TS: {results['heavy_chain']['gdt_ts']:.2f}%")
                print(f"LDDT: {results['heavy_chain']['lddt']:.2f}%")
            
            if light_chain_id and "rmsd" in results["light_chain"]:
                print(f"\nLight Chain ({light_chain_id}) Metrics:")
                print(f"RMSD: {results['light_chain']['rmsd']:.3f} Å")
                print(f"GDT_TS: {results['light_chain']['gdt_ts']:.2f}%")
                print(f"LDDT: {results['light_chain']['lddt']:.2f}%")
                
            if ag_chain_id and "rmsd" in results["antigen_chain"]:
                print(f"\nAntigen Chain ({ag_chain_id}) Metrics:")
                print(f"RMSD: {results['antigen_chain']['rmsd']:.3f} Å")
                print(f"GDT_TS: {results['antigen_chain']['gdt_ts']:.2f}%")
                print(f"LDDT: {results['antigen_chain']['lddt']:.2f}%")
        
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

        fig, axes = plt.subplots(3, 1, figsize=(20, 15), sharey=False)
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
    
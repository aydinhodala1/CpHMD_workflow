def lambda_group(lambda_name: str, group_name: str, atom_charges: list[list[float]], pka: float, dvdl_coeffs: list[float], index: int) -> str:
    '''
    Create the md.mdp CpHMD string for a lambda group

    Inputs:
    lambda_name (str): Name for lambda group in CpHMD
    group_name (str): Name for index group in index.ndx
    atom_charges (2D list of floats): 1st list corresponds to charges at lambda = 0, 2nd corresponds to charges at lambda = 1. Note 1st number in each list must correspond to first atom in index group, second for second atom etc.
    pka (float): pka of calibration state
    dvdl_coeffs (list of floats): Coefficients for polynomial fit of V_corr
    index (int): Index number for lambda group

    Outputs:
    lambda_str (str): Lines for appending to md.mdp file for CpHMD runs. Note assumes all groups are coupled to the same group of buffer particles (charge-restraint-group-index = 1)
    '''
    #Create each line of lambda-dynamics settings
    name_str = f'lambda-dynamics-group-type{index}-name                          = {lambda_name}'
    states_str = f'lambda-dynamics-group-type{index}-n-states                      = 1'
    lam0_str = f'lambda-dynamics-group-type{index}-state-0-charges               =  {str(atom_charges[0])[1:-1].replace(",","")}'
    lam1_str = f'lambda-dynamics-group-type{index}-state-1-charges               =  {str(atom_charges[1])[1:-1].replace(",","")}'
    pka_str = f'lambda-dynamics-group-type{index}-state-1-reference-pka         = {pka}'
    dvdl_str = f'lambda-dynamics-group-type{index}-state-1-dvdl-coefficients     =  {str(dvdl_coeffs)[1:-1].replace(",","")}'
    atom_str = f'lambda-dynamics-atom-set{index}-name                         = {lambda_name}'
    group_str = f'lambda-dynamics-atom-set{index}-index-group-name             = {group_name}'
    barrier_str = f'lambda-dynamics-atom-set{index}-barrier                      = 5.0'
    init_str = f'lambda-dynamics-atom-set{index}-initial-lambda               = 0.5'
    charge_str = f'lambda-dynamics-atom-set{index}-charge-restraint-group-index = 1'
    
    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str, 
                '',atom_str, group_str, barrier_str, init_str, charge_str]
            )

    return lambda_str

print(lambda_group("OLA","lam1",[[3,4,5],[4,5,6]],5,[6],7))

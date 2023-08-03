from opentrons import protocol_api

# this file serves as a template only, replace-from-excel.py and the template excel file is used to change the wells and volumes

metadata = {
    'protocolName': '20230802-130918-ont-protocol.py',
    'author': 'BCL <angel.angelov@kaust.edu.sa>',
    'description': 'ONT plasmid sequencing - normalise templates, add rapid adapter, pool',
    'apiLevel': '2.8'
}

sourcewells1=['A1','B1','C1','D1','E1','F1','G1','H1','A2','B2','C2','D2','E2','F2','G2','H2','A3','B3','C3','D3','E3','F3',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
destwells1=['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3', 'A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4', 'A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5', 'A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6', 'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7', 'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8', 'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9', 'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10', 'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11', 'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12']
volume1=[3.99071231514369, 1.43896100128071, 1.31049051426148, 4.79411793638479, 1.25896261994499, 1.25525452366443, 5, 3.57591176853526, 1.69886529896907, 3.25037577109194, 2.57216650711065, 3.06174898431049, 5, 1.22288193271768, 2.45783640909091, 2.16289604, 1.7030677480315, 2.00268151851852, 2.39788917960089, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#sourcewells2=['A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
destwells2=['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3', 'A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4', 'A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5', 'A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6', 'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7', 'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8', 'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9', 'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10', 'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11', 'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12']
volume2=[1.00928768485631, 3.56103899871929, 3.68950948573852, 0.205882063615206, 3.74103738005501, 3.74474547633557, 0, 1.42408823146474, 3.30113470103093, 1.74962422890806, 2.42783349288935, 1.93825101568951, 0, 3.77711806728232, 2.54216359090909, 2.83710396, 3.2969322519685, 2.99731848148148, 2.60211082039911, 4.5, 4.5, 4.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
sourcewells3=['A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2','G2','H2','A3','B3','C3','D3','E3','F3','G3','H3',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
destwells3=['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3', 'A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4', 'A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5', 'A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6', 'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7', 'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8', 'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9', 'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10', 'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11', 'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12']
volume3=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

watersource = 'A1'
finaltube = 'B1'
# use 1 ul barcode and 11 ul total rxn vol if it is gDNA, for plasmid use half volumes
# the decision is based on the barcode volumes (volume3) 
barcode_vol = 0.5
total_rxn_vol = 5.5

######################## Calculations for full column transfer  - for rapid barcode plate #######################
# the requirement is that:
# for ONE source column all rows go to ONE dest column AND 
# there has to be a row correspondence A-A, B-B...H-H 
# AND the volumes are the same for the whole column

scols3_fulltransfer = []
dcols3_fulltransfer = []

for i in range(0, 95, 8):
    scols3 = [col[1:] for col in sourcewells3[i:i + 8]]
    dcols3 = [col[1:] for col in destwells3[i:i + 8]]

    if([row[:1] for row in sourcewells3[i:i + 8]] ==  [row[:1] for row in destwells3[i:i + 8]] and # there is row correspondence
        scols3.count(scols3[0]) == len(scols3)): # all wells in the batch of 8 are the same column
    # collect data for transfer
        scols3_fulltransfer.append( scols3[0] )
        dcols3_fulltransfer.append( dcols3[0] )

    # set the vol3 for the whole col transfers to 0
    for i, v in enumerate(destwells3):
        if v[1:] in dcols3_fulltransfer:
            volume3[i] = 0

######################## Calculations for full column transfer  - for rapid barcode plate #######################

# exit early if there is something wrong with the dest wells
if len(destwells1) != 96:
    exit("Please make sure that there are 96 destination wells! Check the excel template is correct...")

def run(ctx: protocol_api.ProtocolContext):
    ctx.comment("Starting ONT plasmid sequencing protocol")

    destplate = ctx.load_labware('pcrplate_96_wellplate_200ul', '5', 'Destination plate') # stack of 96 well base plate and PCR plate
    sourceplate = ctx.load_labware('pcrplate_96_wellplate_200ul', '4', 'Source plate') # stack of 96 well base plate and PCR plate
    barcodeplate = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '1', 'Rapid barcode plate')
    sourcetube = ctx.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '7', 'Tube rack')

    tips20_single = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['10', '11']]
    tips20_multi = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['3']]

    s20 = ctx.load_instrument('p20_single_gen2', mount='left', tip_racks=tips20_single)
    m20 = ctx.load_instrument('p20_multi_gen2', mount='right', tip_racks=tips20_multi)

    # distribute water without tip change first
    ctx.comment("================= Starting water transfer ==========================")
    s20.distribute(	
        volume2,
        sourcetube.wells_by_name()[watersource], 
        [ destplate.wells_by_name()[i] for i in destwells2 ], 
        new_tip = 'always', 
        touch_tip = False
    )
    

    ctx.comment("================= Starting plasmid transfer ==========================")
    # add plasmid, changing tip
    
    s20.transfer(
        [v for v in volume1 if v > 0],
        [ sourceplate.wells_by_name()[v] for i, v in enumerate(sourcewells1) if volume1[i] > 0],
        [ destplate.wells_by_name()[v] for i, v in enumerate(destwells1) if volume1[i] > 0], 
        new_tip = 'always',
        mix_after = (1, 3), 
        blow_out = True, 
        blowout_location = 'destination well'
        )
            

    # add barcodes, full columns if possible, has to be as fast as possible
    ctx.comment("================= Starting barcode transfer ==========================")
    
    # pause - this is optional in the Shiny app so that the protocol can be used to just do DNA adjustment to given conc
# optional pause #    ctx.pause("Optional pause before barcode addition") 

    for i, v in enumerate(scols3_fulltransfer):
        ctx.comment("Full column transfer barcode plate: " + str(barcode_vol) + "ul from A" + v + " to A" + dcols3_fulltransfer[i])
        m20.transfer(
        barcode_vol, 
        barcodeplate.wells_by_name()['A' + scols3_fulltransfer[i]], 
        destplate.wells_by_name()['A' + dcols3_fulltransfer[i]], 
        new_tip = 'always', 
        mix_after = (3, 4), 
        blow_out = True, 
        blowout_location = 'destination well' 
        )
        ctx.comment("--------------------------------------")

    # s20 transfer barcodes for what is left
    for i, v in enumerate(volume3):
        if v > 0:
            ctx.comment("s20 transfer barcode plate")
            s20.transfer(
                v,
                barcodeplate.wells_by_name()[sourcewells3[i]], 
                destplate.wells_by_name()[destwells3[i]], 
                new_tip = 'always', 
                mix_after = (3, 4), 
                blow_out = True, 
        		blowout_location = 'destination well'
            )
            ctx.comment("--------------------------------------")
    
    # indicate it is ready
    for _ in range(3):
        ctx.set_rail_lights(False)
        ctx.delay(1)
        ctx.set_rail_lights(True)
        ctx.delay(1)

    ctx.comment("Please incubate the destination plate at 30°C for 2 minutes and 80°C for 2 minutes")
    ctx.comment("When incubations are ready place back the plate on 5")
    ctx.pause("Please incubate the destination plate at 30°C for 2 minutes and 80°C for 2 minutes. When incubations are ready place back the plate on 5")

    # Pool
    ctx.comment("================= Pool samples =========================")

    # use same tip as everything is aspirated

    s20.consolidate(
        total_rxn_vol,
        [ destplate.wells_by_name()[v] for i, v in enumerate(destwells1) if volume1[i] > 0], 
        sourcetube.wells_by_name()[finaltube]
    )

    poolvol = len([v for v in volume1 if v > 0]) * total_rxn_vol

    ctx.comment("Final volume of pool: " + str(poolvol) + " ul")
    ctx.comment("Add " + str(poolvol) + " ul SPRI beads, incubate 5 minutes. Wash beads 2x with 80 EtOH, aspirate all EtOH and dry 30 sec. Resuspend in 15 ul EB, incubate 10 minutes, take out supernatant and continue with RAP addition")
    
    ctx.comment("================= End =========================")
    

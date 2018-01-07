import puzzlotope.Chem as Chem

Elements = (
            Chem.Element('H', (
                [1.0078250322 , 0.9999],
                [2.01410177781, 0.0001],
            )),
            Chem.Element('O', (
                [15.994914620, 0.9976],
                [16.999131757, 0.0004],
                [17.999159613, 0.0020]
            )),
            Chem.Element('C', (
                [12.0        , 0.9894],
                [13.003354835, 0.0106]
            )),
            Chem.Element('Ni', (
                [57.935342, 0.680769],
                [59.930786, 0.262231],
                [60.931056, 0.011399],
                [61.928345, 0.036345],
            )),
            Chem.Element('Mg', (
                [23.98504170, 0.79],
                [24.9858370 , 0.10],
                [25.9825930 , 0.11],
            )),
            Chem.Element('Cl', (
                [34.9688527, 0.76],
                [36.9659026, 0.24]
            )),
        )

Blocks = (
            Chem.Block('Ni'     , [0,1,2,3,4]   ),
#            Chem.Block('Mg'     , [2,]          ),
            Chem.Block('C6H5'   , [-1,]         , tags=[Chem.BlockTag.optional, Chem.BlockTag.canOxydize]),
            Chem.Block('C5H8'   , [0,]          , tags=[Chem.BlockTag.optional]),
#            Chem.Block('Cl'     , [-1,]         , tags=[Chem.BlockTag.optional]),
            Chem.Block('OH'     , [-1,]         , tags=[Chem.BlockTag.optional]),
            Chem.Block('O'      , [-2,]        , tags=[Chem.BlockTag.optional]),
#            Chem.Block('H2O' ,[0,])
        )

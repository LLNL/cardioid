/* Remember that down in the cell models the units don't necessarily
 *  correspond to the internal units of Cardioid.  The units in this map
 *  are the units the cell model expects the variables to have.
*/
HandleMap& OHaraRudy::getHandleMap()
{
  static HandleMap handleMap;
   if (handleMap.size() == 0)
   {
      handleMap["Vm"            ] = CheckpointVarInfo(Vm              , false, "mV"    );
      handleMap["Nai"           ] = CheckpointVarInfo(Nai             , true , "mM"    );
      handleMap["Nass"          ] = CheckpointVarInfo(Nass            , true , "mM"    );
      handleMap["Ki"            ] = CheckpointVarInfo(Ki              , true , "mM"    );
      handleMap["Kss"           ] = CheckpointVarInfo(Kss             , true , "mM"    );
      handleMap["Cai"           ] = CheckpointVarInfo(Cai             , true , "mM"    );
      handleMap["Cass"          ] = CheckpointVarInfo(Cass            , true , "mM"    );
      handleMap["Cansr"         ] = CheckpointVarInfo(Cansr           , true , "mM"    );
      handleMap["Cajsr"         ] = CheckpointVarInfo(Cajsr           , true , "mM"    );
      handleMap["m"             ] = CheckpointVarInfo(m               , true , "1"     );
      handleMap["hFast"         ] = CheckpointVarInfo(hFast           , true , "1"     );
      handleMap["hSlow"         ] = CheckpointVarInfo(hSlow           , true , "1"     );
      handleMap["j"             ] = CheckpointVarInfo(j               , true , "1"     );
      handleMap["hCaMKSlow"     ] = CheckpointVarInfo(hCaMKSlow       , true , "1"     );
      handleMap["jCaMK"         ] = CheckpointVarInfo(jCaMK           , true , "1"     );
      handleMap["mL"            ] = CheckpointVarInfo(mL              , true , "1"     );
      handleMap["hL"            ] = CheckpointVarInfo(hL              , true , "1"     );
      handleMap["hLCaMK"        ] = CheckpointVarInfo(hLCaMK          , true , "1"     );
      handleMap["a"             ] = CheckpointVarInfo(a               , true , "1"     );
      handleMap["iFast"         ] = CheckpointVarInfo(iFast           , true , "1"     );
      handleMap["iSlow"         ] = CheckpointVarInfo(iSlow           , true , "1"     );
      handleMap["aCaMK"         ] = CheckpointVarInfo(aCaMK           , true , "1"     );
      handleMap["iCaMKFast"     ] = CheckpointVarInfo(iCaMKFast       , true , "1"     );
      handleMap["iCaMKSlow"     ] = CheckpointVarInfo(iCaMKSlow       , true , "1"     );
      handleMap["d"             ] = CheckpointVarInfo(d               , true , "1"     );
      handleMap["fFast"         ] = CheckpointVarInfo(fFast           , true , "1"     );
      handleMap["fSlow"         ] = CheckpointVarInfo(fSlow           , true , "1"     );
      handleMap["fCaFast"       ] = CheckpointVarInfo(fCaFast         , true , "1"     );
      handleMap["fCaSlow"       ] = CheckpointVarInfo(fCaSlow         , true , "1"     );
      handleMap["jCa"           ] = CheckpointVarInfo(jCa             , true , "1"     );
      handleMap["n"             ] = CheckpointVarInfo(n               , true , "1"     );
      handleMap["fCaMKFast"     ] = CheckpointVarInfo(fCaMKFast       , true , "1"     );
      handleMap["fCaCaMKFast"   ] = CheckpointVarInfo(fCaCaMKFast     , true , "1"     );
      handleMap["XrFast"        ] = CheckpointVarInfo(XrFast          , true , "1"     );
      handleMap["XrSlow"        ] = CheckpointVarInfo(XrSlow          , true , "1"     );
      handleMap["Xs1"           ] = CheckpointVarInfo(Xs1             , true , "1"     );
      handleMap["Xs2"           ] = CheckpointVarInfo(Xs2             , true , "1"     );
      handleMap["XK1"           ] = CheckpointVarInfo(XK1             , true , "1"     );
      handleMap["JrelNP"        ] = CheckpointVarInfo(JrelNP          , true , "mM/ms" );
      handleMap["JrelCaMK"      ] = CheckpointVarInfo(JrelCaMK        , true , "mM/ms" );
      handleMap["CaMKtrap"      ] = CheckpointVarInfo(CaMKtrap        , true , "1"     );
      handleMap["GNaFast"       ] = CheckpointVarInfo(GNaFast         , false, "mS/uF" );
      handleMap["GNaL"          ] = CheckpointVarInfo(GNaL            , false, "mS/uF" );
      handleMap["Gto"           ] = CheckpointVarInfo(Gto             , false, "mS/uF" );
      handleMap["PCaL"          ] = CheckpointVarInfo(PCaL            , false, "cm/s"  );
      handleMap["PCaNa"         ] = CheckpointVarInfo(PCaNa           , false, "cm/s"  );
      handleMap["PCaK"          ] = CheckpointVarInfo(PCaK            , false, "cm/s"  );
      handleMap["GKr"           ] = CheckpointVarInfo(GKr             , false, "mS/uF" );
      handleMap["GK1"           ] = CheckpointVarInfo(GK1             , false, "mS/uF" );
      handleMap["GKs"           ] = CheckpointVarInfo(GKs             , false, "mS/uF" );
      handleMap["GNaCai"        ] = CheckpointVarInfo(GNaCai          , false, "uA/uF" );
      handleMap["GNaCass"       ] = CheckpointVarInfo(GNaCass         , false, "uA/uF" );
      handleMap["PNaK"          ] = CheckpointVarInfo(PNaK            , false, "mV/mM" );
      handleMap["PNab"          ] = CheckpointVarInfo(PNab            , false, "cm/s"  );
      handleMap["PCab"          ] = CheckpointVarInfo(PCab            , false, "cm/s"  );
      handleMap["GKb"           ] = CheckpointVarInfo(GKb             , false, "mS/uF" );
      handleMap["GpCa"          ] = CheckpointVarInfo(GpCa            , false, "mS/uF" );
      handleMap["alphaJrelNP"   ] = CheckpointVarInfo(alphaJrelNP     , false, "mM/mV" );
      handleMap["betaJrelNP"    ] = CheckpointVarInfo(betaJrelNP      , false, "mM/mS" );
      handleMap["alphaJrelCaMK" ] = CheckpointVarInfo(alphaJrelCaMK   , false, "mM/mV" );
      handleMap["betaJrelCaMK"  ] = CheckpointVarInfo(betaJrelCaMK    , false, "mM/mS" );
      handleMap["cJup"          ] = CheckpointVarInfo(cJup            , false, "mM/ms" );
      handleMap["CMDN"          ] = CheckpointVarInfo(CMDN            , false, "mM"    );
      handleMap["aDelta"        ] = CheckpointVarInfo(aDelta          , false, "1"     );
      assert(handleMap.size() == nVars);
   }
   return handleMap;
}

void OHaraRudy::setValue(int varHandle, double value)
{
 
   switch (varHandle)
   {
      case undefinedName   : assert(false)                   ; break;
      case Vm              : state_.Vm               = value; break;
      case Nai             : state_.Nai              = value; break;
      case Nass            : state_.Nass             = value; break;
      case Ki              : state_.Ki               = value; break;
      case Kss             : state_.Kss              = value; break;
      case Cai             : state_.Cai              = value; break;
      case Cass            : state_.Cass             = value; break;
      case Cansr           : state_.Cansr            = value; break;
      case Cajsr           : state_.Cajsr            = value; break;
      case m               : state_.m                = value; break;
      case hFast           : state_.hFast            = value; break;
      case hSlow           : state_.hSlow            = value; break;
      case j               : state_.j                = value; break;
      case hCaMKSlow       : state_.hCaMKSlow        = value; break;
      case jCaMK           : state_.jCaMK            = value; break;
      case mL              : state_.mL               = value; break;
      case hL              : state_.hL               = value; break;
      case hLCaMK          : state_.hLCaMK           = value; break;
      case a               : state_.a                = value; break;
      case iFast           : state_.iFast            = value; break;
      case iSlow           : state_.iSlow            = value; break;
      case aCaMK           : state_.aCaMK            = value; break;
      case iCaMKFast       : state_.iCaMKFast        = value; break;
      case iCaMKSlow       : state_.iCaMKSlow        = value; break;
      case d               : state_.d                = value; break;
      case fFast           : state_.fFast            = value; break;
      case fSlow           : state_.fSlow            = value; break;
      case fCaFast         : state_.fCaFast          = value; break;
      case fCaSlow         : state_.fCaSlow          = value; break;
      case jCa             : state_.jCa              = value; break;
      case n               : state_.n                = value; break;
      case fCaMKFast       : state_.fCaMKFast        = value; break;
      case fCaCaMKFast     : state_.fCaCaMKFast      = value; break;
      case XrFast          : state_.XrFast           = value; break;
      case XrSlow          : state_.XrSlow           = value; break;
      case Xs1             : state_.Xs1              = value; break;
      case Xs2             : state_.Xs2              = value; break;
      case XK1             : state_.XK1              = value; break;
      case JrelNP          : state_.JrelNP           = value; break;
      case JrelCaMK        : state_.JrelCaMK         = value; break;
      case CaMKtrap        : state_.CaMKtrap         = value; break;
      case GNaFast         : cellParms_->GNaFast          = value; break;
      case GNaL            : cellParms_->GNaL             = value; break;
      case Gto             : cellParms_->Gto              = value; break;
      case PCaL            : cellParms_->PCaL             = value; break;
      case PCaNa           : cellParms_->PCaNa            = value; break;
      case PCaK            : cellParms_->PCaK             = value; break;
      case GKr             : cellParms_->GKr              = value; break;
      case GK1             : cellParms_->GK1              = value; break;
      case GKs             : cellParms_->GKs              = value; break;
      case GNaCai          : cellParms_->GNaCai           = value; break;
      case GNaCass         : cellParms_->GNaCass          = value; break;
      case PNaK            : cellParms_->PNaK             = value; break;
      case PNab            : cellParms_->PNab             = value; break;
      case PCab            : cellParms_->PCab             = value; break;
      case GKb             : cellParms_->GKb              = value; break;
      case GpCa            : cellParms_->GpCa             = value; break;
      case alphaJrelNP     : cellParms_->alphaJrelNP      = value; break;
      case betaJrelNP      : cellParms_->betaJrelNP       = value; break;
      case alphaJrelCaMK   : cellParms_->alphaJrelCaMK    = value; break;
      case betaJrelCaMK    : cellParms_->betaJrelCaMK     = value; break;
      case cJup            : cellParms_->cJup             = value; break;
      case CMDN            : cellParms_->CMDN             = value; break;
      case aDelta          : cellParms_->aDelta           = value; break;
      case nVars           : assert(false)                  ; break;
   }
}

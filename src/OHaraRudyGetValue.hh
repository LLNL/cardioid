double OHaraRudy::getValue(int varHandle) const
{
 
   switch (varHandle)
   {
      case undefinedName   : assert(false)                 ; break;
      case Vm              : return state_.Vm              ; break;
      case Nai             : return state_.Nai             ; break;
      case Nass            : return state_.Nass            ; break;
      case Ki              : return state_.Ki              ; break;
      case Kss             : return state_.Kss             ; break;
      case Cai             : return state_.Cai             ; break;
      case Cass            : return state_.Cass            ; break;
      case Cansr           : return state_.Cansr           ; break;
      case Cajsr           : return state_.Cajsr           ; break;
      case m               : return state_.m               ; break;
      case hFast           : return state_.hFast           ; break;
      case hSlow           : return state_.hSlow           ; break;
      case j               : return state_.j               ; break;
      case hCaMKSlow       : return state_.hCaMKSlow       ; break;
      case jCaMK           : return state_.jCaMK           ; break;
      case mL              : return state_.mL              ; break;
      case hL              : return state_.hL              ; break;
      case hLCaMK          : return state_.hLCaMK          ; break;
      case a               : return state_.a               ; break;
      case iFast           : return state_.iFast           ; break;
      case iSlow           : return state_.iSlow           ; break;
      case aCaMK           : return state_.aCaMK           ; break;
      case iCaMKFast       : return state_.iCaMKFast       ; break;
      case iCaMKSlow       : return state_.iCaMKSlow       ; break;
      case d               : return state_.d               ; break;
      case fFast           : return state_.fFast           ; break;
      case fSlow           : return state_.fSlow           ; break;
      case fCaFast         : return state_.fCaFast         ; break;
      case fCaSlow         : return state_.fCaSlow         ; break;
      case jCa             : return state_.jCa             ; break;
      case n               : return state_.n               ; break;
      case fCaMKFast       : return state_.fCaMKFast       ; break;
      case fCaCaMKFast     : return state_.fCaCaMKFast     ; break;
      case XrFast          : return state_.XrFast          ; break;
      case XrSlow          : return state_.XrSlow          ; break;
      case Xs1             : return state_.Xs1             ; break;
      case Xs2             : return state_.Xs2             ; break;
      case XK1             : return state_.XK1             ; break;
      case JrelNP          : return state_.JrelNP          ; break;
      case JrelCaMK        : return state_.JrelCaMK        ; break;
      case CaMKtrap        : return state_.CaMKtrap        ; break;
      case GNaFast         : return cellParms_->GNaFast         ; break;
      case GNaL            : return cellParms_->GNaL            ; break;
      case Gto             : return cellParms_->Gto             ; break;
      case PCaL            : return cellParms_->PCaL            ; break;
      case PCaNa           : return cellParms_->PCaNa           ; break;
      case PCaK            : return cellParms_->PCaK            ; break;
      case GKr             : return cellParms_->GKr             ; break;
      case GK1             : return cellParms_->GK1             ; break;
      case GKs             : return cellParms_->GKs             ; break;
      case GNaCai          : return cellParms_->GNaCai          ; break;
      case GNaCass         : return cellParms_->GNaCass         ; break;
      case PNaK            : return cellParms_->PNaK            ; break;
      case PNab            : return cellParms_->PNab            ; break;
      case PCab            : return cellParms_->PCab            ; break;
      case GKb             : return cellParms_->GKb             ; break;
      case GpCa            : return cellParms_->GpCa            ; break;
      case alphaJrelNP     : return cellParms_->alphaJrelNP     ; break;
      case betaJrelNP      : return cellParms_->betaJrelNP      ; break;
      case alphaJrelCaMK   : return cellParms_->alphaJrelCaMK   ; break;
      case betaJrelCaMK    : return cellParms_->betaJrelCaMK    ; break;
      case cJup            : return cellParms_->cJup            ; break;
      case CMDN            : return cellParms_->CMDN            ; break;
      case aDelta          : return cellParms_->aDelta          ; break;
      case nVars           : assert(false)                  ; break;
   }
   return 0.;
}

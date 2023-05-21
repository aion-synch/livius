#ifndef CORENC_METHOD_DG_FLUX_H_
#define CORENC_METHOD_DG_FLUX_H_
namespace corenc
{
	namespace method
	{
		enum class DGFlux
		{
			EIP,
			EBaumannOden,
			EBaumannOdenIP,
			ENIPG,
			EUpwind,
			ECentral,
			ELaxFriedrichs,
			IIP,
			IBaumannOden,
			IBaumannOdenIP,
			INIPG,
			IUpwind,
			ICentral,
			ILaxFriedrichs,
			CUSTOM,
			NOFLUX,
		};

	}
}
#endif // !CORENC_METHOD_DG_FLUX_H_

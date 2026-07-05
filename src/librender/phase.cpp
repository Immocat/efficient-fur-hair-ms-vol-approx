#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "PhaseFunctionSamplingRecord[" << endl
		<< "  mRec = " << indent(mRec.toString()) << "," << endl
		<< "  wi = " << wi.toString() << "," << endl
		<< "  wo = " << wo.toString() << "," << endl
		<< "  mode = " << mode << endl
		<< "]";
	return oss.str();
}

void PhaseFunction::configure() {
	m_type = 0;
}

Float PhaseFunction::pdf(const PhaseFunctionSamplingRecord &pRec) const {
	return eval(pRec);
}

bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return false;
}

Float PhaseFunction::sigmaDir(Float cosTheta) const {
	Log(EError, "%s::sigmaDir(Float) is not implemented (this is not "
				"an anisotropic medium!)",
		getClass()->getName().c_str());
	return 0.0f;
}

Float PhaseFunction::sigmaDirMax() const {
	Log(EError, "%s::sigmaDirMax() is not implemented (this is not "
				"an anisotropic medium!)",
		getClass()->getName().c_str());
	return 0.0f;
}

Float PhaseFunction::getMeanCosine() const {
	Log(EError, "%s::getMeanCosine() is not implemented!",
		getClass()->getName().c_str());
	return 0.0f;
}

Spectrum PhaseFunction::albedoDir(const Vector &d, const Vector &orientation) const {
	Log(EError, "%s::albedoDir(Float) is not implemented (this is not "
				"an anisotropic medium!)",
		getClass()->getName().c_str());
	return Spectrum(0.0f);
}

bool PhaseFunction::needsDirectionallyVaryingAlbedo() const {
	return false;
}

Float PhaseFunction::sigmaDir(const Vector &d, const Vector &orientation) const {
	return sigmaDir(dot(d, orientation));
}

DynamicMatrix PhaseFunction::PrecomputeHairMatrix(int bands, int nSamples) const {
	Log(EError, "%s::PrecomputeHairMatrix() is not implemented",
		getClass()->getName().c_str());
}

MatrixRGB PhaseFunction::PrecomputeHairMatrixRGB(int bands, int nSamples) const {
	Log(EError, "%s::PrecomputeHairMatrix() is not implemented",
		getClass()->getName().c_str());
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(SpectrumPhaseFunction, true, PhaseFunction)
MTS_NAMESPACE_END

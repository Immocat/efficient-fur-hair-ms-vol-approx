/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/shvector.h>

MTS_NAMESPACE_BEGIN

class TestSphericalHarmonics : public TestCase {
public:
    MTS_BEGIN_TESTCASE()
    MTS_DECLARE_TEST(test01_shRotation)
    MTS_DECLARE_TEST(test02_shRotationIdentity)
    MTS_DECLARE_TEST(test03_shRotationOrthonormality)
    MTS_DECLARE_TEST(test04_shRotationComposition)
    MTS_DECLARE_TEST(test05_shRotationFirstOrderBasis)
    MTS_DECLARE_TEST(test06_shFirstOrderConvention)
    MTS_DECLARE_TEST(test07_shEvalHelpersConsistency)
    MTS_DECLARE_TEST(test08_shVectorRGBProject)
    MTS_DECLARE_TEST(test09_rotateSHRGBIdentity)
    MTS_DECLARE_TEST(test10_rotateSHRGBMatchesPerChannelRotation)
    MTS_DECLARE_TEST(test11_shSampler)
    MTS_END_TESTCASE()

    void test01_shRotation() {
        /* Generate a random SH expansion, rotate it and
           spot-check 100 times against the original evaluated
           at appropriately rotated positions */

        ref<Random> random = new Random();
        int bands = 8;

        SHVector vec1(bands);
        for (int l=0; l<bands; ++l)
            for (int m=-l; m<=l; ++m)
                vec1(l, m) = random->nextFloat();

        Vector axis(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat())));
        Transform trafo = Transform::rotate(axis, random->nextFloat()*360);
        SHRotation rot(vec1.getBands());

        SHVector::rotation(trafo, rot);
        SHVector vec2(bands);

        rot(vec1, vec2);

        for (int i=0; i<100; ++i) {
            Vector dir1(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat()))), dir2;
            trafo(dir1, dir2);

            Float value1 = vec1.eval(dir2);
            Float value2 = vec2.eval(dir1);
            assertEqualsEpsilon(value1, value2, Epsilon);
        }
    }

    void test02_shRotationIdentity() {
        const int bands = 6;
        SHRotation rot(bands);
        SHVector::rotation(Transform(), rot);

        for (int l = 0; l < bands; ++l) {
            const SHRotation::Matrix &M = rot.blocks[l];
            for (int i = 0; i < M.rows(); ++i) {
                for (int j = 0; j < M.cols(); ++j) {
                    const Float expected = (i == j) ? (Float) 1 : (Float) 0;
                    assertEqualsEpsilon(M(i, j), expected, 1e-5f);
                }
            }
        }
    }

    void test03_shRotationOrthonormality() {
        ref<Random> random = new Random();
        const int bands = 6;

        for (int iter = 0; iter < 10; ++iter) {
            Vector axis(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat())));
            Transform trafo = Transform::rotate(axis, random->nextFloat() * 360);
            SHRotation rot(bands);
            SHVector::rotation(trafo, rot);

            for (int l = 0; l < bands; ++l) {
                const SHRotation::Matrix &M = rot.blocks[l];
                SHRotation::Matrix shouldBeIdentity = M * M.transpose();
                for (int i = 0; i < shouldBeIdentity.rows(); ++i) {
                    for (int j = 0; j < shouldBeIdentity.cols(); ++j) {
                        const Float expected = (i == j) ? (Float) 1 : (Float) 0;
                        assertEqualsEpsilon(shouldBeIdentity(i, j), expected, 1e-4f);
                    }
                }
            }
        }
    }

    void test04_shRotationComposition() {
        ref<Random> random = new Random();
        const int bands = 5;

        for (int iter = 0; iter < 10; ++iter) {
            Vector axis1(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat())));
            Vector axis2(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat())));
            Transform t1 = Transform::rotate(axis1, random->nextFloat() * 360);
            Transform t2 = Transform::rotate(axis2, random->nextFloat() * 360);
            Transform combined = t2 * t1;

            SHRotation rot1(bands), rot2(bands), rotCombined(bands);
            SHVector::rotation(t1, rot1);
            SHVector::rotation(t2, rot2);
            SHVector::rotation(combined, rotCombined);

            for (int l = 0; l < bands; ++l) {
                SHRotation::Matrix product = rot1.blocks[l] * rot2.blocks[l];
                const SHRotation::Matrix &expected = rotCombined.blocks[l];
                for (int i = 0; i < expected.rows(); ++i) {
                    for (int j = 0; j < expected.cols(); ++j) {
                        assertEqualsEpsilon(product(i, j), expected(i, j), 2e-4f);
                    }
                }
            }
        }
    }

    void test05_shRotationFirstOrderBasis() {
        const int bands = 2;
        const Vector directions[] = {
            Vector(1, 0, 0),
            Vector(0, 1, 0),
            Vector(0, 0, 1),
            normalize(Vector(1, 2, 3))
        };
        const Transform rotations[] = {
            Transform::rotate(Vector(0, 0, 1), 90.0f),
            Transform::rotate(Vector(0, 1, 0), 90.0f),
            Transform::rotate(Vector(1, 0, 0), 90.0f)
        };

        for (size_t rIdx = 0; rIdx < sizeof(rotations) / sizeof(rotations[0]); ++rIdx) {
            SHRotation rot(bands);
            SHVector::rotation(rotations[rIdx], rot);

            Log(EInfo, "test05: rotation[%i] first-order block:", (int) rIdx);
            for (int row = 0; row < rot.blocks[1].rows(); ++row) {
                Log(EInfo, "  row %i = [%f, %f, %f]",
                    row,
                    rot.blocks[1](row, 0),
                    rot.blocks[1](row, 1),
                    rot.blocks[1](row, 2));
            }

            for (int basisIdx = 0; basisIdx < 3; ++basisIdx) {
                SHVector source(bands), target(bands);
                source.clear();
                source(1, basisIdx - 1) = 1.0f;
                rot(source, target);

                for (size_t dIdx = 0; dIdx < sizeof(directions) / sizeof(directions[0]); ++dIdx) {
                    Vector rotatedDir;
                    rotations[rIdx](directions[dIdx], rotatedDir);
                    Float expected = source.eval(rotatedDir);
                    Float actual = target.eval(directions[dIdx]);
                    if (std::abs(actual - expected) > 1e-5f) {
                        Log(EInfo,
                            "test05 mismatch: rotation=%i basis(l=1,m=%i) dir=(%f,%f,%f) rotatedDir=(%f,%f,%f) expected=%f actual=%f",
                            (int) rIdx, basisIdx - 1,
                            directions[dIdx].x, directions[dIdx].y, directions[dIdx].z,
                            rotatedDir.x, rotatedDir.y, rotatedDir.z,
                            expected, actual);
                    }
                    assertEqualsEpsilon(actual, expected, 1e-5f);
                }
            }
        }
    }

    void test06_shFirstOrderConvention() {
        const int bands = 2;
        const Vector directions[] = {
            Vector(1, 0, 0),
            Vector(0, 1, 0),
            Vector(0, 0, 1),
            normalize(Vector(1, 2, 3))
        };
        const char *labels[] = { "+X", "+Y", "+Z", "normalize(1,2,3)" };

        for (int basisIdx = 0; basisIdx < 3; ++basisIdx) {
            SHVector basis(bands);
            basis.clear();
            basis(1, basisIdx - 1) = 1.0f;

            Log(EInfo, "test06: basis coefficient (l=1,m=%i)", basisIdx - 1);
            for (size_t dIdx = 0; dIdx < sizeof(directions) / sizeof(directions[0]); ++dIdx) {
                Log(EInfo, "  %s -> %f", labels[dIdx], basis.eval(directions[dIdx]));
            }
        }
    }

    void test07_shEvalHelpersConsistency() {
        const int bands = 4;
        const int coeffCount = bands * bands;
        const Vector directions[] = {
            Vector(1, 0, 0),
            Vector(0, 1, 0),
            Vector(0, 0, 1),
            normalize(Vector(1, 2, 3)),
            normalize(Vector(-2, 1, 4))
        };

        for (size_t dIdx = 0; dIdx < sizeof(directions) / sizeof(directions[0]); ++dIdx) {
            std::vector<Float> stdValues;
            DynamicVector eigenValues;
            EvalSH_All(bands, directions[dIdx], stdValues);
            EvalSH_All_Eigen(bands, directions[dIdx], eigenValues);

            for (int index = 0; index < coeffCount; ++index) {
                const Float singleValue = EvalSH(index, directions[dIdx]);
                assertEqualsEpsilon(singleValue, stdValues[index], 1e-6f);
                assertEqualsEpsilon(singleValue, eigenValues[index], 1e-6f);
            }
        }
    }

    void test08_shVectorRGBProject() {
        const int bands = 4;
        const Vector dir = normalize(Vector(1, 2, 3));
        Spectrum energy;
        energy.fromLinearRGB(0.25f, 0.5f, 0.75f);

        SHVectorRGB rgb(bands);
        SHVector basisValues(bands);
        rgb.project(dir, energy, basisValues);

        for (int l = 0; l < bands; ++l) {
            for (int m = -l; m <= l; ++m) {
                const int index = l * (l + 1) + m;
                const Float basis = EvalSH(index, dir);
                assertEqualsEpsilon(basisValues(l, m), basis, 1e-6f);
                assertEqualsEpsilon(rgb.r(l, m), energy[0] * basis, 1e-6f);
                assertEqualsEpsilon(rgb.g(l, m), energy[1] * basis, 1e-6f);
                assertEqualsEpsilon(rgb.b(l, m), energy[2] * basis, 1e-6f);
            }
        }
    }

    void test09_rotateSHRGBIdentity() {
        const int bands = 4;
        SHVectorRGB source(bands), target;

        for (int l = 0; l < bands; ++l) {
            for (int m = -l; m <= l; ++m) {
                const Float base = (Float) (l * 10 + m);
                source.r(l, m) = base + 0.1f;
                source.g(l, m) = base + 0.2f;
                source.b(l, m) = base + 0.3f;
            }
        }

        SHRotation rot(bands);
        SHVector::rotation(Transform(), rot);
        RotateSHRGB(rot, source, target);

        for (int l = 0; l < bands; ++l) {
            for (int m = -l; m <= l; ++m) {
                assertEqualsEpsilon(target.r(l, m), source.r(l, m), 1e-6f);
                assertEqualsEpsilon(target.g(l, m), source.g(l, m), 1e-6f);
                assertEqualsEpsilon(target.b(l, m), source.b(l, m), 1e-6f);
            }
        }
    }

    void test10_rotateSHRGBMatchesPerChannelRotation() {
        const int bands = 4;
        const Vector axis = normalize(Vector(1, -2, 3));
        const Transform trafo = Transform::rotate(axis, 37.0f);

        SHVectorRGB source(bands), target;
        for (int l = 0; l < bands; ++l) {
            for (int m = -l; m <= l; ++m) {
                const Float base = (Float) (l * 10 + m);
                source.r(l, m) = base + 1.0f;
                source.g(l, m) = base + 2.0f;
                source.b(l, m) = base + 3.0f;
            }
        }

        SHRotation rot(bands);
        SHVector::rotation(trafo, rot);
        RotateSHRGB(rot, source, target);

        SHVector expectedR(bands), expectedG(bands), expectedB(bands);
        rot(source.r, expectedR);
        rot(source.g, expectedG);
        rot(source.b, expectedB);

        for (int l = 0; l < bands; ++l) {
            for (int m = -l; m <= l; ++m) {
                assertEqualsEpsilon(target.r(l, m), expectedR(l, m), 1e-6f);
                assertEqualsEpsilon(target.g(l, m), expectedG(l, m), 1e-6f);
                assertEqualsEpsilon(target.b(l, m), expectedB(l, m), 1e-6f);
            }
        }
    }

    struct ClampedCos {
        Vector axis;
        ClampedCos(Vector axis) : axis(axis) { }
        Float operator()(const Vector &w) const { return std::max((Float) 0, dot(w, axis)); }
    };

    void test11_shSampler() {
        /* Draw 100 samples from a SH expansion of a clamped cosine-shaped
           distribution and verify the returned probabilities */
        int bands = 13, numSamples = 100, depth = 12;

        Vector v = normalize(Vector(1, 2, 3));
        ref<Random> random = new Random();
        SHVector clampedCos = SHVector(bands);
        clampedCos.project(ClampedCos(v), numSamples);
        //Float clampedCosError = clampedCos.l2Error(ClampedCos(v), numSamples);
        clampedCos.normalize();

        //cout << "Projection error = " << clampedCosError << endl;
        //cout << "Precomputing mip-maps" << endl;
        ref<SHSampler> sampler = new SHSampler(bands, depth);
        //cout << "Done: "<< sampler->toString() << endl;
        Float accum = 0;
        int nsamples = 100, nInAvg = 0;
        for (int i=0; i<=nsamples; ++i) {
            Point2 sample(random->nextFloat(), random->nextFloat());
            Float pdf1 = sampler->warp(clampedCos, sample);
            Float pdf2 = dot(v, sphericalDirection(sample.x, sample.y))/M_PI;
            Float relerr = std::abs(pdf1-pdf2)/pdf2;
            if (pdf2 > 0.01) {
                accum += relerr; ++nInAvg;
                assertTrue(relerr < 0.08);
            }
        }
        assertTrue(accum / nInAvg < 0.01);
    }
};

MTS_EXPORT_TESTCASE(TestSphericalHarmonics, "Testcase for Spherical Harmonics code")
MTS_NAMESPACE_END

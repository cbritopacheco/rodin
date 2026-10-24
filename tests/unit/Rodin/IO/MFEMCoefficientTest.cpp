/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Independent MFEM coefficient layouts and native input diagnostics.
#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/IO.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  template <class Scalar>
  class MFEMCoefficients : public ::testing::Test
  {
    protected:
      /// @brief Gives distinct exactly representable entries for every node and component.
      static Scalar entry(size_t node, size_t component)
      {
        const Real real = 1 + 8 * Real(node) + Real(component);
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(real, -2 - 4 * Real(node) - Real(component));
        else
          return real;
      }

      /// @brief Builds an independent MFEM stream in either supported component ordering.
      static std::string payload(
        size_t nodes, size_t components, size_t ordering, bool discontinuous)
      {
        std::ostringstream stream;
        stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
               << "# independent coefficients\nFiniteElementSpace\n"
               << "FiniteElementCollection: " << (discontinuous ? "L2_1D_P0" : "H1_1D_P1")
               << "\nVDim: " << components << "\nOrdering: " << ordering << "\n\n";
        for (size_t outer = 0; outer < (ordering == 0 ? components : nodes); ++outer)
        {
          for (size_t inner = 0; inner < (ordering == 0 ? nodes : components); ++inner)
          {
            stream << entry(ordering == 0 ? inner : outer, ordering == 0 ? outer : inner)
                   << '\n';
          }
        }
        return stream.str();
      }

      /// @brief Checks independent layouts and printer output, including rectangular matrix components.
      template <size_t Rank, bool Discontinuous>
      void verify()
      {
        using Range = std::conditional_t<Rank == 0, Scalar,
          std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
            Math::SpatialMatrix<Scalar>>>;
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
        using FES = std::conditional_t<Discontinuous, P0<Range>, P1<Range>>;
        auto space = [&] {
          if constexpr (Rank == 0)
            return FES(mesh);
          else if constexpr (Rank == 1)
            return FES(mesh, 3);
          else
            return FES(mesh, 2, 3);
        }();
        const size_t nodes = Discontinuous ? mesh.getCellCount() : mesh.getVertexCount();
        const size_t components = space.getVectorDimension();
        GridFunction field(space);
        GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Scalar>> loader(field);
        for (size_t ordering : {size_t(0), size_t(1)})
        {
          SCOPED_TRACE(ordering);
          field.getData().setConstant(Scalar(-99));
          std::istringstream input(payload(nodes, components, ordering, Discontinuous));
          ASSERT_NO_THROW(loader.load(input));
          for (size_t node = 0; node < nodes; ++node)
          {
            for (size_t component = 0; component < components; ++component)
            {
              const size_t index = (Rank == 2 || Discontinuous)
                ? node * components + component
                : node + component * nodes;
              EXPECT_EQ(field.getData()[index], entry(node, component));
            }
          }
          std::ostringstream output;
          output << std::setprecision(std::numeric_limits<Real>::max_digits10);
          GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Scalar>>(field).print(
            output);
          std::istringstream printed(output.str());
          std::string line;
          for (size_t i = 0; i < 4; ++i)
          {
            ASSERT_TRUE(bool(std::getline(printed, line)));
            if (i == 2)
            {
              EXPECT_EQ(line, "VDim: " + std::to_string(components));
            }
            if (i == 3)
            {
              EXPECT_EQ(line, "Ordering: 0");
            }
          }
          Scalar actual;
          for (size_t component = 0; component < components; ++component)
          {
            for (size_t node = 0; node < nodes; ++node)
            {
              ASSERT_TRUE(bool(printed >> actual));
              EXPECT_EQ(actual, entry(node, component));
            }
          }
          EXPECT_FALSE(bool(printed >> actual));
        }
      }

      /// @brief Wrong P0 coefficient counts preserve scalar/vector data and permit reuse.
      template <size_t Rank>
      void verifyMalformedP0Coefficients()
      {
        using Range = std::conditional_t<Rank == 0, Scalar, Math::SpatialVector<Scalar>>;
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
        auto space = [&] {
          if constexpr (Rank == 0)
            return P0<Range>(mesh);
          else
            return P0<Range>(mesh, 3);
        }();
        GridFunction field(space);
        GridFunctionLoader<FileFormat::MFEM, decltype(space), Math::Vector<Scalar>>
          loader(field);
        const size_t nodes = mesh.getCellCount();
        const size_t components = space.getVectorDimension();
        const size_t count = nodes * components;
        for (const size_t ordering : {size_t(0), size_t(1)})
        {
          SCOPED_TRACE(ordering);
          const std::string valid = payload(nodes, components, ordering, true);
          const std::string header = valid.substr(0, valid.find("\n\n") + 2);
          std::vector<std::pair<std::string, size_t>> invalid{{header + "bad\n", 0}};
          std::ostringstream partial;
          partial << header << entry(0, 0) << "\nbad\n";
          invalid.emplace_back(partial.str(), 1);
          for (const size_t found : {size_t(0), count - 1, count + 1})
          {
            std::ostringstream stream;
            stream << header;
            for (size_t i = 0; i < found; ++i)
              stream << entry(i % nodes, i / nodes) << '\n';
            invalid.emplace_back(stream.str(), found);
          }
          for (const auto& [data, found] : invalid)
          {
            SCOPED_TRACE(data);
            for (size_t i = 0; i < count; ++i)
              field.getData()[i] = -entry(i % nodes, i / nodes);
            const auto before = field.getData().eval();
            std::istringstream input(data);
            try
            {
              loader.load(input);
              FAIL() << "Wrong P0 coefficient count was accepted";
            }
            catch (const Alert::Exception& error)
            {
              const std::string message = error.what();
              EXPECT_NE(message.find("In member function"), std::string::npos);
              EXPECT_NE(message.find("load"), std::string::npos);
              EXPECT_NE(message.find("Expected " + std::to_string(count) +
                          " coefficient values, but found " + std::to_string(found)),
                std::string::npos);
            }
            EXPECT_EQ(field.getData(), before);
            for (const size_t recoveryOrdering : {size_t(0), size_t(1)})
            {
              std::istringstream recovery(
                payload(nodes, components, recoveryOrdering, true));
              ASSERT_NO_THROW(loader.load(recovery));
              for (size_t node = 0; node < nodes; ++node)
              {
                for (size_t component = 0; component < components; ++component)
                {
                  EXPECT_EQ(field.getData()[node * components + component],
                    entry(node, component));
                }
              }
            }
          }
        }
      }

      /// @brief Rejects malformed matrix headers and coefficients before modifying data.
      template <bool Discontinuous>
      void verifyMalformedMatrix()
      {
        using Range = Math::SpatialMatrix<Scalar>;
        using FES = std::conditional_t<Discontinuous, P0<Range>, P1<Range>>;
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
        FES space(mesh, 2, 3);
        GridFunction field(space);
        GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Scalar>> loader(field);
        const size_t nodes = space.getScalarSpace().getSize();
        const size_t components = space.getVectorDimension();
        const auto header = [&](
                              const std::string& dimension, const std::string& ordering) {
          return std::string("FiniteElementSpace\nFiniteElementCollection: ") +
            (Discontinuous ? "L2_1D_P0" : "H1_1D_P1") + "\n" + dimension + "\n" +
            ordering + "\n\n";
        };
        std::vector<std::pair<std::string, std::string>> invalid;
        for (const std::string dimension :
          {"Dimension: 6", "VDim: bad", "VDim: 5", "VDim:"})
        {
          invalid.emplace_back(
            header(dimension, "Ordering: 0"), "Invalid MFEM matrix field header");
        }
        for (const std::string ordering :
          {"Order: 0", "Ordering: bad", "Ordering: 2", "Ordering:"})
        {
          invalid.emplace_back(
            header("VDim: 6", ordering), "Invalid MFEM matrix field header");
        }
        for (const size_t ordering : {size_t(0), size_t(1)})
        {
          for (const size_t validValues : {size_t(0), size_t(1), nodes * components - 1})
          {
            std::ostringstream data;
            data << header("VDim: 6", "Ordering: " + std::to_string(ordering));
            for (size_t i = 0; i < validValues; ++i)
              data << Scalar(1 + i) << '\n';
            invalid.emplace_back(data.str(), "Truncated MFEM matrix data");
            invalid.emplace_back(data.str() + "bad\n", "Truncated MFEM matrix data");
          }
        }
        for (const auto& [data, diagnostic] : invalid)
        {
          SCOPED_TRACE(data);
          for (Index i = 0; i < static_cast<size_t>(field.getData().size()); ++i)
          {
            if constexpr (std::is_same_v<Scalar, Complex>)
              field.getData()[i] = Scalar(-17 - i, 3 + i);
            else
              field.getData()[i] = Scalar(-17 - i);
          }
          const auto before = field.getData().eval();
          std::istringstream input(data);
          try
          {
            loader.load(input);
            FAIL() << "Malformed matrix input was accepted";
          }
          catch (const Alert::Exception& error)
          {
            const std::string message = error.what();
            EXPECT_NE(message.find("member function"), std::string::npos);
            EXPECT_NE(message.find("load"), std::string::npos);
            EXPECT_NE(message.find(diagnostic), std::string::npos);
          }
          catch (const std::exception& error)
          {
            ADD_FAILURE() << "Expected a native member-function diagnostic: "
                          << error.what();
          }
          EXPECT_EQ(field.getData(), before);
          for (const size_t ordering : {size_t(0), size_t(1)})
          {
            std::istringstream valid(payload(nodes, components, ordering, Discontinuous));
            ASSERT_NO_THROW(loader.load(valid));
            for (size_t node = 0; node < nodes; ++node)
            {
              for (size_t component = 0; component < components; ++component)
              {
                EXPECT_EQ(
                  field.getData()[node * components + component], entry(node, component));
              }
            }
          }
        }
      }
  };
  using MFEMScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(MFEMCoefficients, MFEMScalars);

  /// @brief Default output streams preserve nonintegral real and complex coefficients.
  TYPED_TEST(MFEMCoefficients, DefaultStreamPreservesCoefficientPrecision)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    const auto verify = [&]<class FES>(const FES& space) {
      GridFunction<FES, Math::Vector<TypeParam>> field(space);
      for (size_t i = 0; i < space.getSize(); ++i)
      {
        const Real real = Real(i + 1) / 7;
        if constexpr (std::is_same_v<TypeParam, Complex>)
          field.getData()[i] = TypeParam(real, -Real(i + 2) / 11);
        else
          field.getData()[i] = real;
      }
      std::ostringstream output;
      GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<TypeParam>>(field).print(
        output);
      GridFunction<FES, Math::Vector<TypeParam>> restored(space);
      restored.getData().setConstant(TypeParam(-99));
      std::istringstream input(output.str());
      GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<TypeParam>>(restored).load(
        input);
      EXPECT_EQ(restored.getData(), field.getData());
      input >> std::ws;
      EXPECT_TRUE(input.eof());
    };
    verify(P0<TypeParam>(mesh));
    verify(P1<TypeParam>(mesh));
    verify(H1<3, TypeParam>(std::integral_constant<size_t, 3>{}, mesh));
  }

  TYPED_TEST(MFEMCoefficients, P1ScalarBothOrderings)
  {
    this->template verify<0, false>();
  }
  TYPED_TEST(MFEMCoefficients, P1VectorBothOrderings)
  {
    this->template verify<1, false>();
  }
  TYPED_TEST(MFEMCoefficients, P1MatrixBothOrderings)
  {
    this->template verify<2, false>();
  }
  TYPED_TEST(MFEMCoefficients, P0ScalarBothOrderings)
  {
    this->template verify<0, true>();
  }
  TYPED_TEST(MFEMCoefficients, P0VectorBothOrderings)
  {
    this->template verify<1, true>();
  }
  TYPED_TEST(MFEMCoefficients, P0MatrixBothOrderings)
  {
    this->template verify<2, true>();
  }

  TYPED_TEST(MFEMCoefficients, P1MalformedMatrixPreservesDestinationAndRecovers)
  {
    this->template verifyMalformedMatrix<false>();
  }
  TYPED_TEST(MFEMCoefficients, P0MalformedMatrixPreservesDestinationAndRecovers)
  {
    this->template verifyMalformedMatrix<true>();
  }

  TYPED_TEST(MFEMCoefficients, P0MalformedScalarAndVectorPreserveDestinationAndRecover)
  {
    this->template verifyMalformedP0Coefficients<0>();
    this->template verifyMalformedP0Coefficients<1>();
  }

  /// @brief Invalid first and truncated later P1 coefficients raise native diagnostics before changing the field.
  TYPED_TEST(MFEMCoefficients, P1MalformedCoefficientsPreserveDestination)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    P1<TypeParam> space(mesh);
    GridFunction field(space);
    GridFunctionLoader<FileFormat::MFEM, decltype(space), Math::Vector<TypeParam>> loader(
      field);
    for (const std::string data : {"bad\n2\n3\n4\n", "1\n2\n", "1\nbad\n3\n4\n"})
    {
      SCOPED_TRACE(data);
      field.getData().setConstant(TypeParam(-17));
      const auto before = field.getData().eval();
      std::istringstream input("FiniteElementSpace\nFiniteElementCollection: H1_1D_P1\n"
                               "VDim: 1\nOrdering: 0\n\n" +
        data);
      try
      {
        loader.load(input);
        FAIL() << "Malformed coefficient input was accepted";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("member function"), std::string::npos);
        EXPECT_NE(message.find("load"), std::string::npos);
      }
      catch (const std::exception& error)
      {
        ADD_FAILURE() << "Expected a native member-function diagnostic: " << error.what();
      }
      EXPECT_EQ(field.getData(), before);
    }
  }
}

/// @brief Header validation and truncated H1 input preserve existing coefficient storage.
TEST(MFEMFieldValidation, HeaderAndPayload)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  const auto check = [&]<class Space>(Space& space, const std::string& collection) {
    using Field = GridFunction<Space, Math::Vector<Real>>;
    Field field(space);
    GridFunctionLoader<FileFormat::MFEM, Space, Math::Vector<Real>> loader(field);
    const std::string header =
      "FiniteElementSpace\nFiniteElementCollection: " + collection +
      "\nVDim: 1\nOrdering: 0\n\n";
    std::ostringstream coefficients;
    for (size_t i = 0; i < space.getSize(); ++i)
      coefficients << i + 1 << ' ';
    const std::string valid = header + coefficients.str();
    std::istringstream packed(valid);
    ASSERT_NO_THROW(loader.load(packed));
    const auto before = field.getData().eval();
    const auto reject = [&](std::string input) {
      std::istringstream stream(input);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      EXPECT_EQ(field.getData(), before);
    };
    reject("BadSpace\n" + valid.substr(valid.find('\n') + 1));
    for (const auto& [from, to] :
      std::vector<std::pair<std::string, std::string>>{{"VDim: 1", "VDim: 0"},
        {"VDim: 1", "VDim: 2"}, {"Ordering: 0", "Ordering: 2"}, {collection, "H1_3D_P9"}})
    {
      std::string invalid = valid;
      invalid.replace(invalid.find(from), from.size(), to);
      reject(invalid);
    }
    reject(header + "1");
    std::istringstream recovery(valid);
    EXPECT_NO_THROW(loader.load(recovery));
  };
  P1<Real> linear(mesh);
  H1<2, Real> quadratic(std::integral_constant<size_t, 2>{}, mesh);
  check(linear, "H1_1D_P1");
  check(quadratic, "H1_1D_P2");
}

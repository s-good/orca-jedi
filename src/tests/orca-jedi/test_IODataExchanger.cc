/*
 * (C) British Crown Copyright 2023 Met Office
 */

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "atlas/parallel/mpi/mpi.h"

#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace.h"

#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "orca-jedi/geometry/IODataExchanger.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

  CASE ("Test IODataExchanger") {

    atlas::OrcaGrid grid("ORCA2_T");
    auto meshgen_config = grid.meshgenerator();

    atlas::MeshGenerator meshgen(meshgen_config);
    auto partitioner_config = grid.partitioner();
    partitioner_config.set("type", "checkerboard");
    auto partitioner = atlas::grid::Partitioner(partitioner_config);
    auto mesh = meshgen.generate(grid, partitioner);
    auto funcSpace = atlas::functionspace::NodeColumns(mesh);

    int nlevels = 2;
    IODataExchanger exchanger(grid, nlevels);

    const atlas::OrcaGrid orcaGrid = atlas::OrcaGrid(mesh.grid());
    int nx = orcaGrid.nx() + orcaGrid.haloEast() + orcaGrid.haloWest();
    int ny = orcaGrid.ny() + orcaGrid.haloNorth() + orcaGrid.haloSouth();

    // setup the initial field
    std::cout << "setup the initial field" << std::endl;
    atlas::Field field(funcSpace.createField<double>(
                          atlas::option::name("testfield") |
                          atlas::option::levels(nlevels)));
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());

    std::cout << " field_view.shape(0) " << field_view.shape(0) << std::endl;
    for (size_t k =0; k < nlevels; ++k) {
      for (size_t i = 0; i < field_view.shape(0); ++i) {
        if (ghost(i)) continue;
        field_view(i, k) = 0;
      }
    }
    field.haloExchange();

    // create a global field
    std::cout << "create a global field" << std::endl;
    atlas::Field globalField = exchanger.gather(field);

    // fill a buffer with unique identifiers by creating a serially distributed mesh
    std::vector<double> buffer;
    {
      auto global_meshgen_config = grid.meshgenerator();
      atlas::MeshGenerator global_meshgen(global_meshgen_config);
      auto global_partitioner_config = grid.partitioner();
      global_partitioner_config.set("type", "serial");
      auto global_partitioner = atlas::grid::Partitioner(global_partitioner_config);
      auto global_mesh = global_meshgen.generate(grid, partitioner);
      auto master_global_index = atlas::array::make_view<atlas::gidx_t, 1>(global_mesh.nodes().field("master_global_index"));
      auto global_ij = atlas::array::make_view<int32_t, 2>(global_mesh.nodes().field("ij"));
      buffer.resize(nlevels*nx*ny);
      for (size_t k = 0; k < nlevels; ++k) {
        for (size_t iNode = 0; iNode < master_global_index.shape(0); ++k) {
            const int i = global_ij(iNode, 0) + orcaGrid.haloWest();
            const int j = global_ij(iNode, 1) + orcaGrid.haloSouth();
            buffer[k*nx*ny + j*ny + i] = master_global_index(iNode);
        }
      }
    }

    // copy that buffer into a global_view
    std::cout << "copy that buffer into a global_view" << std::endl;
    if (atlas::mpi::comm().rank() == 0) {
      auto global_view = atlas::array::make_view<double, 2>(globalField);
      std::cout << " global_view.shape(0) " << global_view.shape(0)
                << " nx " << nx
                << " ny " << ny
                << std::endl;
      for (size_t k = 0; k < nlevels; ++k) {
        for (size_t iNode = 0; iNode < global_view.shape(0); ++iNode) {
          global_view(iNode, k) = buffer[exchanger.globalToBufferIndex(iNode)];
        }
      }
    }

    // clear the buffer contents
    buffer.clear();

    // scatter the global field
    std::cout << "scatter the global field" << std::endl;
    exchanger.scatter(globalField, field);

    field.haloExchange();

    auto local_master_global_index = atlas::array::make_view<atlas::gidx_t, 1>(mesh.nodes().field("master_global_index"));
    // check the results
    for (size_t k = 0; k < nlevels; ++k) {
      for (size_t iNode = 0; iNode < field_view.shape(0); ++iNode) {
        EXPECT_EQUAL(field_view(iNode, k), local_master_global_index(iNode));
      }
    }
  }

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}

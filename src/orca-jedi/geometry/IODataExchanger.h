#pragma once

#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/field.h"
#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/array/ArrayView.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

class ORCAGlobalNodeToBufferIndex {
  private:
    int32_t ix_glb_max;
    int32_t iy_glb_max;
    int32_t glbarray_offset;
    int32_t glbarray_jstride;
    int32_t nx_halo_WE;
    int32_t ny_halo_NS;
  public:
    std::vector<int32_t> mapping;

    ORCAGlobalNodeToBufferIndex() {};
    explicit ORCAGlobalNodeToBufferIndex(const atlas::OrcaGrid& orcaGrid,
                                         const atlas::array::ArrayView<int32_t, 2> ij) {
        // calculate the offset and stride for conversion between i,j
        // coordinates and a flat file buffer
        iy_glb_max = orcaGrid.ny() + orcaGrid.haloNorth() - 1;
        ix_glb_max = orcaGrid.nx() + orcaGrid.haloEast() - 1;
        nx_halo_WE = orcaGrid.nx() + orcaGrid.haloEast() + orcaGrid.haloWest();
        ny_halo_NS = orcaGrid.ny() + orcaGrid.haloNorth()
          + orcaGrid.haloSouth();
        int iy_glb_min = -orcaGrid.haloSouth();
        int ix_glb_min = -orcaGrid.haloWest();
        glbarray_offset  = -(nx_halo_WE * iy_glb_min) - ix_glb_min;
        glbarray_jstride = nx_halo_WE;

        // create a mapping between the global node index and the buffer file index
        for (int32_t iNode=0; iNode < ij.shape(0); ++iNode) {
          int32_t i = ij(iNode, 0);
          int32_t j = ij(iNode, 1);
          ATLAS_ASSERT(i <= ix_glb_max,
              std::to_string(i) + " > " + std::to_string(ix_glb_max));
          ATLAS_ASSERT(j <= iy_glb_max,
              std::to_string(j) + " > " + std::to_string(iy_glb_max));
          mapping.emplace_back(glbarray_offset + j * glbarray_jstride + i);
        }
    }

};

    //struct IndexGlbArray {
    //    int32_t ix_glb_max;
    //    int32_t iy_glb_max;
    //    int32_t glbarray_offset;
    //    int32_t glbarray_jstride;
    //    int32_t nx_halo_WE;
    //    int32_t ny_halo_NS;

    //    explicit IndexGlbArray(const atlas::OrcaGrid& orcaGrid) {
    //        iy_glb_max = orcaGrid.ny() + orcaGrid.haloNorth() - 1;
    //        ix_glb_max = orcaGrid.nx() + orcaGrid.haloEast() - 1;

    //        nx_halo_WE = orcaGrid.nx() + orcaGrid.haloEast() + orcaGrid.haloWest();
    //        ny_halo_NS = orcaGrid.ny() + orcaGrid.haloNorth()
    //          + orcaGrid.haloSouth();

    //        // vector of local indices: necessary for remote indices of ghost nodes
    //        int iy_glb_min = -orcaGrid.haloSouth();
    //        int ix_glb_min = -orcaGrid.haloWest();
    //        glbarray_offset  = -(nx_halo_WE * iy_glb_min) - ix_glb_min;
    //        glbarray_jstride = nx_halo_WE;
    //    }

    //    int32_t operator()(int32_t i, int32_t j) {
    //        ATLAS_ASSERT(i <= ix_glb_max,
    //            std::to_string(i) + " > " + std::to_string(ix_glb_max));
    //        ATLAS_ASSERT(j <= iy_glb_max,
    //            std::to_string(j) + " > " + std::to_string(iy_glb_max));
    //        return glbarray_offset + j * glbarray_jstride + i;
    //    }
    //};

    //std::vector<size_t> request_indices_;
    //std::vector<std::vector<size_t>> remote_indices_;
    //atlas::array::ArrayView<int32_t, 1> ghost;

    //int32_t nx_;
    //int32_t ny_;
    //int32_t nlevels_;

    //explicit IODataExchanger(const atlas::Mesh& mesh, const int32_t nlevels) :
    //    ghost{atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost())}
    //{
    //    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
    //    const atlas::OrcaGrid orcaGrid = atlas::OrcaGrid(mesh.grid());
    //    IndexGlbArray index_glbarray(orcaGrid);
    //    nx_ = index_glbarray.nx_halo_WE;
    //    ny_ = index_glbarray.ny_halo_NS;
    //    nlevels_ = nlevels;
    //    //const size_t numNodes = mesh.nodes().size();
    //    //atlas_omp_for(int k = 0; k < nlevels_; ++k) {
    //    //  for (size_t inode = 0; inode < numNodes; ++inode) {
    //    //    if (ghost(inode)) continue;
    //    //    request_indices_.emplace_back(k*nx_*ny_ + index_glbarray(ij(inode, 0), ij(inode, 1)));
    //    //  }
    //    //}

    //    // send/receive remote_indices_ data required from each rank
    //    //remote_indices_.resize( atlas::mpi::size() );
    //    //atlas::mpi::comm().allToAll(request_indices_, remote_indices_);
    //}

    //void scatter(const std::vector<double>& buffer,
    //             atlas::array::ArrayView<double, 2>& field_view) {

    //    std::vector<double> local_buffer(request_indices_.size());

    //    // this isn't scatter, but it does need to send out the data according
    //    // to remote indices. what should this be? some kind of intelligent
    //    // scatter?
    //    //atlas::mpi::comm().scatter(buffer, local_buffer, 0);

    //    const size_t numNodes = field_view.shape(0);
    //    size_t j = 0;
    //    atlas_omp_for(int k = 0; k < nlevels_; ++k) {
    //      for (size_t inode = 0; inode < numNodes; ++inode) {
    //        if (ghost(inode)) continue;
    //        field_view(inode, k) = local_buffer[j];
    //        ++j;
    //      }
    //    }
    //}
class IODataExchanger {

 public:
    explicit IODataExchanger(const atlas::Grid& grid, const int32_t n_levels) :
      grid_(grid), n_levels_(n_levels)
    {
    //  // create a global mesh just to extract the ij field for long enough to
    //  // produce a mapping from the ORCA ij mesh to the file buffer indices
    //  // TODO: Future optimisation this can probably be done more efficiently
    //  // using OrcaGrid::index2ij, however my understanding is that this function
    //  // doesn't return exactly the same thing at the moment.
    //  auto meshgen_config = grid.meshgenerator();
    //  atlas::MeshGenerator meshgen(meshgen_config);
    //  auto partitioner_config = grid.partitioner();
    //  partitioner_config.set("type", "serial");
    //  auto partitioner = atlas::grid::Partitioner(partitioner_config);
    //  auto mesh = meshgen.generate(grid, partitioner);
    //  auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    //  const atlas::OrcaGrid orcaGrid = atlas::OrcaGrid(grid);
    //  globalToBufferMap = ORCAGlobalNodeToBufferIndex(orcaGrid, ij);
    }

    void scatter(atlas::Field& globalField, atlas::Field& localField) {
    //  auto& functionSpace = globalField.functionspace();
    //  functionSpace.scatter(globalField, localField);
    }

    atlas::Field gather(const atlas::Field& field) {
    //  if (field.metadata().get<bool>("global") == false) {
    //    atlas::array::DataType atlasType = field.datatype();
    //    atlas::util::Config atlasOptions = atlas::option::name(field.name()) |
    //                                       atlas::option::levels(field.levels()) |
    //                                       atlas::option::datatype(atlasType) |
    //                                       atlas::option::global(0);
    //    const auto& functionSpace = field.functionspace();
    //    atlas::Field globalField = functionSpace.createField(atlasOptions);
    //    field.haloExchange();
    //    std::cout << "functionSpace.gather(field, globalField)" << std::endl;
    //    functionSpace.gather(field, globalField);
    //    return globalField;
    //  } else {
        return field;
    //  }
    }

    int32_t globalToBufferIndex (int32_t iNode) const { return -1; } // globalToBufferMap.mapping[iNode]; }

  private:
    ORCAGlobalNodeToBufferIndex globalToBufferMap;
    int32_t n_levels_;
    atlas::Grid grid_;
};
}  // namespace orcamodel

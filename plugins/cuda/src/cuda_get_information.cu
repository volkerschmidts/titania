
#include <cuda_get_information.hpp>
#include <iostream>


void cuda_print_information()
{
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties (&deviceProp, 0);

   int i;

   std::cout << "ECCEnabled  : " << deviceProp.ECCEnabled << std::endl;
   std::cout << "asyncEngineCount  : " << deviceProp.asyncEngineCount << std::endl;
   std::cout << "canMapHostMemory  : " << deviceProp.canMapHostMemory << std::endl;
   std::cout << "canUseHostPointerForRegisteredMem  : " << deviceProp.canUseHostPointerForRegisteredMem << std::endl;
   std::cout << "clockRate  : " << deviceProp.clockRate << std::endl;
   std::cout << "computeMode  : " << deviceProp.computeMode << std::endl;
   std::cout << "computePreemptionSupported  : " << deviceProp.computePreemptionSupported << std::endl;
   std::cout << "concurrentKernels  : " << deviceProp.concurrentKernels << std::endl;
   std::cout << "concurrentManagedAccess  : " << deviceProp.concurrentManagedAccess << std::endl;
   std::cout << "cooperativeLaunch  : " << deviceProp.cooperativeLaunch << std::endl;
   std::cout << "cooperativeMultiDeviceLaunch  : " << deviceProp.cooperativeMultiDeviceLaunch << std::endl;
   std::cout << "deviceOverlap  : " << deviceProp.deviceOverlap << std::endl;
   std::cout << "directManagedMemAccessFromHost  : " << deviceProp.directManagedMemAccessFromHost << std::endl;
   std::cout << "globalL1CacheSupported  : " << deviceProp.globalL1CacheSupported << std::endl;
   std::cout << "hostNativeAtomicSupported  : " << deviceProp.hostNativeAtomicSupported << std::endl;
   std::cout << "integrated  : " << deviceProp.integrated << std::endl;
   std::cout << "isMultiGpuBoard  : " << deviceProp.isMultiGpuBoard << std::endl;
   std::cout << "kernelExecTimeoutEnabled  : " << deviceProp.kernelExecTimeoutEnabled << std::endl;
   std::cout << "l2CacheSize  : " << deviceProp.l2CacheSize / ( 1024.0)<< std::endl;
   std::cout << "localL1CacheSupported  : " << deviceProp.localL1CacheSupported << std::endl;
   std::cout << "luid  : " << deviceProp.luid << std::endl;
   std::cout << "luidDeviceNodeMask  : " << deviceProp.luidDeviceNodeMask << std::endl;
   std::cout << "major  : " << deviceProp.major << std::endl;
   std::cout << "managedMemory  : " << deviceProp.managedMemory << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxGridSize  : " << deviceProp.maxGridSize[i] << std::endl;

   std::cout << "maxSurface1D  : " << deviceProp.maxSurface1D << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxSurface1DLayered  : " << deviceProp.maxSurface1DLayered[i] << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxSurface2D  : " << deviceProp.maxSurface2D[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxSurface2DLayered  : " << deviceProp.maxSurface2DLayered[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxSurface3D  : " << deviceProp.maxSurface3D[i] << std::endl;

   std::cout << "maxSurfaceCubemap  : " << deviceProp.maxSurfaceCubemap << std::endl;

   for ( i = 0; i < i; ++i )
   std::cout << "maxSurfaceCubemapLayered  : " << deviceProp.maxSurfaceCubemapLayered[i] << std::endl;

   std::cout << "maxTexture1D  : " << deviceProp.maxTexture1D << std::endl;

   for ( i = 0; i < i; ++i )
   std::cout << "maxTexture1DLayered  : " << deviceProp.maxTexture1DLayered[i] << std::endl;

   std::cout << "maxTexture1DLinear  : " << deviceProp.maxTexture1DLinear << std::endl;
   std::cout << "maxTexture1DMipmap  : " << deviceProp.maxTexture1DMipmap << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxTexture2D  : " << deviceProp.maxTexture2D[i] << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxTexture2DGather  : " << deviceProp.maxTexture2DGather[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxTexture2DLayered  : " << deviceProp.maxTexture2DLayered[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxTexture2DLinear  : " << deviceProp.maxTexture2DLinear[i] << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxTexture2DMipmap  : " << deviceProp.maxTexture2DMipmap[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxTexture3D  : " << deviceProp.maxTexture3D[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxTexture3DAlt  : " << deviceProp.maxTexture3DAlt[i] << std::endl;

   std::cout << "maxTextureCubemap  : " << deviceProp.maxTextureCubemap << std::endl;

   for ( i = 0; i < 2; ++i )
   std::cout << "maxTextureCubemapLayered  : " << deviceProp.maxTextureCubemapLayered[i] << std::endl;

   for ( i = 0; i < 3; ++i )
   std::cout << "maxThreadsDim  : " << deviceProp.maxThreadsDim[i] << std::endl;

   std::cout << "maxThreadsPerBlock  : " << deviceProp.maxThreadsPerBlock << std::endl;
   std::cout << "maxThreadsPerMultiProcessor  : " << deviceProp.maxThreadsPerMultiProcessor << std::endl;

   std::cout << "memPitch  : " << deviceProp.memPitch << std::endl;
   std::cout << "memoryBusWidth  : " << deviceProp.memoryBusWidth << std::endl;
   std::cout << "memoryClockRate  : " << deviceProp.memoryClockRate << std::endl;
   std::cout << "minor  : " << deviceProp.minor << std::endl;
   std::cout << "multiGpuBoardGroupID  : " << deviceProp.multiGpuBoardGroupID << std::endl;
   std::cout << "multiProcessorCount  : " << deviceProp.multiProcessorCount << std::endl;
   std::cout << "name  : " << deviceProp.name << std::endl;
   std::cout << "pageableMemoryAccess  : " << deviceProp.pageableMemoryAccess << std::endl;
   std::cout << "pageableMemoryAccessUsesHostPageTables  : " << deviceProp.pageableMemoryAccessUsesHostPageTables << std::endl;
   std::cout << "pciBusID  : " << deviceProp.pciBusID << std::endl;
   std::cout << "pciDeviceID  : " << deviceProp.pciDeviceID << std::endl;
   std::cout << "pciDomainID  : " << deviceProp.pciDomainID << std::endl;
   std::cout << "regsPerBlock  : " << deviceProp.regsPerBlock << std::endl;
   std::cout << "regsPerMultiprocessor  : " << deviceProp.regsPerMultiprocessor << std::endl;
   std::cout << "sharedMemPerBlock  : " << deviceProp.sharedMemPerBlock << std::endl;
   std::cout << "sharedMemPerBlockOptin  : " << deviceProp.sharedMemPerBlockOptin << std::endl;
   std::cout << "sharedMemPerMultiprocessor  : " << deviceProp.sharedMemPerMultiprocessor << std::endl;
   std::cout << "singleToDoublePrecisionPerfRatio  : " << deviceProp.singleToDoublePrecisionPerfRatio << std::endl;
   std::cout << "streamPrioritiesSupported  : " << deviceProp.streamPrioritiesSupported << std::endl;
   std::cout << "surfaceAlignment  : " << deviceProp.surfaceAlignment << std::endl;
   std::cout << "tccDriver  : " << deviceProp.tccDriver << std::endl;
   std::cout << "textureAlignment  : " << deviceProp.textureAlignment << std::endl;
   std::cout << "texturePitchAlignment  : " << deviceProp.texturePitchAlignment << std::endl;
   std::cout << "totalConstMem  : " << deviceProp.totalConstMem << std::endl;
   std::cout << "totalGlobalMem  : " << deviceProp.totalGlobalMem / ( 1024.0*1024.0 )<< std::endl;
   std::cout << "unifiedAddressing  : " << deviceProp.unifiedAddressing << std::endl;
   //std::cout << "uuid  : " << deviceProp.uuid << std::endl;
   std::cout << "warpSize  : " << deviceProp.warpSize << std::endl;
}


void cuda_print_memory_information()
{
   constexpr double bytes_2_mb = ( 1024.0*1024.0 );
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties (&deviceProp, 0);

   size_t free, total;
   //std::cout << "Device name: " << deviceProp.name << std::endl;
   cudaMemGetInfo(&free,&total);
   std::cout << "Free memory: " << free/ bytes_2_mb << " MB\nTotal memory: " << total/ bytes_2_mb << " MB" << std::endl;
}

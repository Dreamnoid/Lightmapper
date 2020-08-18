#define EMBREE_LIB "embree3.dll"
#include <rtcore_buffer.h>
#include <rtcore.h>

namespace embree
{
	typedef RTCDevice (*NewDeviceFunc)(const char* config);
	typedef void (*ReleaseDeviceFunc)(RTCDevice device);
	typedef void (*SetDeviceErrorFunctionFunc)(RTCDevice device, RTCErrorFunction error, void* userPtr);
	typedef ssize_t (*GetDevicePropertyFunc)(RTCDevice device, enum RTCDeviceProperty prop);
	typedef RTCScene (*NewSceneFunc)(RTCDevice device);
	typedef void (*ReleaseSceneFunc)(RTCScene scene);
	typedef unsigned int (*AttachGeometryFunc)(RTCScene scene, RTCGeometry geometry);
	typedef void (*CommitSceneFunc)(RTCScene scene);
	typedef RTCGeometry (*NewGeometryFunc)(RTCDevice device, enum RTCGeometryType type);
	typedef void (*ReleaseGeometryFunc)(RTCGeometry geometry);
	typedef void (*SetSharedGeometryBufferFunc)(RTCGeometry geometry, enum RTCBufferType type, unsigned int slot, enum RTCFormat format, const void* ptr, size_t byteOffset, size_t byteStride, size_t itemCount);
	typedef void (*CommitGeometryFunc)(RTCGeometry geometry);
	typedef void (*Intersect1Func)(RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit* rayhit);
	typedef void (*Intersect4Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit4* rayhit);
	typedef void (*Intersect8Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit8* rayhit);
	typedef void (*Intersect16Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit16* rayhit);
	typedef void (*Occluded1Func)(RTCScene scene, struct RTCIntersectContext* context, struct RTCRay* ray);
	typedef void (*Occluded4Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay4* ray);
	typedef void (*Occluded8Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay8* ray);
	typedef void (*Occluded16Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay16* ray);
	NewDeviceFunc NewDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDeviceErrorFunctionFunc SetDeviceErrorFunction;
	GetDevicePropertyFunc GetDeviceProperty;
	NewSceneFunc NewScene;
	ReleaseSceneFunc ReleaseScene;
	AttachGeometryFunc AttachGeometry;
	CommitSceneFunc CommitScene;
	NewGeometryFunc NewGeometry;
	ReleaseGeometryFunc ReleaseGeometry;
	SetSharedGeometryBufferFunc SetSharedGeometryBuffer;
	CommitGeometryFunc CommitGeometry;
	Intersect1Func Intersect1;
	Intersect4Func Intersect4;
	Intersect8Func Intersect8;
	Intersect16Func Intersect16;
	Occluded1Func Occluded1;
	Occluded4Func Occluded4;
	Occluded8Func Occluded8;
	Occluded16Func Occluded16;

	void init()
	{
		void* embreeLibrary = loadLibrary("embree3.dll");
		if (!embreeLibrary)
		{
			fatalError("embree not installed");
		}

		embree::NewDevice = (embree::NewDeviceFunc)loadSymbol(embreeLibrary, "rtcNewDevice");
		embree::ReleaseDevice = (embree::ReleaseDeviceFunc)loadSymbol(embreeLibrary, "rtcReleaseDevice");
		embree::SetDeviceErrorFunction = (embree::SetDeviceErrorFunctionFunc)loadSymbol(embreeLibrary, "rtcSetDeviceErrorFunction");
		embree::GetDeviceProperty = (embree::GetDevicePropertyFunc)loadSymbol(embreeLibrary, "rtcGetDeviceProperty");
		embree::NewScene = (embree::NewSceneFunc)loadSymbol(embreeLibrary, "rtcNewScene");
		embree::ReleaseScene = (embree::ReleaseSceneFunc)loadSymbol(embreeLibrary, "rtcReleaseScene");
		embree::AttachGeometry = (embree::AttachGeometryFunc)loadSymbol(embreeLibrary, "rtcAttachGeometry");
		embree::CommitScene = (embree::CommitSceneFunc)loadSymbol(embreeLibrary, "rtcCommitScene");
		embree::NewGeometry = (embree::NewGeometryFunc)loadSymbol(embreeLibrary, "rtcNewGeometry");
		embree::ReleaseGeometry = (embree::ReleaseGeometryFunc)loadSymbol(embreeLibrary, "rtcReleaseGeometry");
		embree::SetSharedGeometryBuffer = (embree::SetSharedGeometryBufferFunc)loadSymbol(embreeLibrary, "rtcSetSharedGeometryBuffer");
		embree::CommitGeometry = (embree::CommitGeometryFunc)loadSymbol(embreeLibrary, "rtcCommitGeometry");
		embree::Intersect1 = (embree::Intersect1Func)loadSymbol(embreeLibrary, "rtcIntersect1");
		embree::Intersect4 = (embree::Intersect4Func)loadSymbol(embreeLibrary, "rtcIntersect4");
		embree::Intersect8 = (embree::Intersect8Func)loadSymbol(embreeLibrary, "rtcIntersect8");
		embree::Intersect16 = (embree::Intersect16Func)loadSymbol(embreeLibrary, "rtcIntersect16");
		embree::Occluded1 = (embree::Occluded1Func)loadSymbol(embreeLibrary, "rtcOccluded1");
		embree::Occluded4 = (embree::Occluded4Func)loadSymbol(embreeLibrary, "rtcOccluded4");
		embree::Occluded8 = (embree::Occluded8Func)loadSymbol(embreeLibrary, "rtcOccluded8");
		embree::Occluded16 = (embree::Occluded16Func)loadSymbol(embreeLibrary, "rtcOccluded16");
	}
};
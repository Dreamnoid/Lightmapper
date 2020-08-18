#define OIDN_LIB "OpenImageDenoise.dll"

namespace oidn
{
	typedef OIDNDevice (*NewDeviceFunc)(OIDNDeviceType type);
	typedef void (*CommitDeviceFunc)(OIDNDevice device);
	typedef void (*ReleaseDeviceFunc)(OIDNDevice device);
	typedef void (*SetDevice1bFunc)(OIDNDevice device, const char* name, bool value);
	typedef void (*SetDeviceErrorFunctionFunc)(OIDNDevice device, OIDNErrorFunction func, void* userPtr);
	typedef OIDNFilter (*NewFilterFunc)(OIDNDevice device, const char* type);
	typedef void (*SetFilterProgressMonitorFunctionFunc)(OIDNFilter filter, OIDNProgressMonitorFunction func, void* userPtr);
	typedef void (*SetSharedFilterImageFunc)(OIDNFilter filter, const char* name, void* ptr, OIDNFormat format, size_t width, size_t height, size_t byteOffset, size_t bytePixelStride, size_t byteRowStride);
	typedef void (*SetFilter1bFunc)(OIDNFilter filter, const char* name, bool value);
	typedef void (*CommitFilterFunc)(OIDNFilter filter);
	typedef void (*ExecuteFilterFunc)(OIDNFilter filter);
	typedef void (*ReleaseFilterFunc)(OIDNFilter filter);
	NewDeviceFunc NewDevice;
	CommitDeviceFunc CommitDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDevice1bFunc SetDevice1b;
	SetDeviceErrorFunctionFunc SetDeviceErrorFunction;
	NewFilterFunc NewFilter;
	SetFilterProgressMonitorFunctionFunc SetFilterProgressMonitorFunction;
	SetSharedFilterImageFunc SetSharedFilterImage;
	SetFilter1bFunc SetFilter1b;
	CommitFilterFunc CommitFilter;
	ExecuteFilterFunc ExecuteFilter;
	ReleaseFilterFunc ReleaseFilter;

	void init()
	{
		void* oidnLibrary = loadLibrary(OIDN_LIB);
		if(!oidnLibrary)
		{
			fatalError("OIDN not installed");
		}

		oidn::NewDevice = (oidn::NewDeviceFunc)loadSymbol(oidnLibrary, "oidnNewDevice");
		oidn::CommitDevice = (oidn::CommitDeviceFunc)loadSymbol(oidnLibrary, "oidnCommitDevice");
		oidn::ReleaseDevice = (oidn::ReleaseDeviceFunc)loadSymbol(oidnLibrary, "oidnReleaseDevice");
		oidn::SetDevice1b = (oidn::SetDevice1bFunc)loadSymbol(oidnLibrary, "oidnSetDevice1b");
		oidn::SetDeviceErrorFunction = (oidn::SetDeviceErrorFunctionFunc)loadSymbol(oidnLibrary, "oidnSetDeviceErrorFunction");
		oidn::NewFilter = (oidn::NewFilterFunc)loadSymbol(oidnLibrary, "oidnNewFilter");
		oidn::SetFilterProgressMonitorFunction = (oidn::SetFilterProgressMonitorFunctionFunc)loadSymbol(oidnLibrary, "oidnSetFilterProgressMonitorFunction");
		oidn::SetSharedFilterImage = (oidn::SetSharedFilterImageFunc)loadSymbol(oidnLibrary, "oidnSetSharedFilterImage");
		oidn::SetFilter1b = (oidn::SetFilter1bFunc)loadSymbol(oidnLibrary, "oidnSetFilter1b");
		oidn::CommitFilter = (oidn::CommitFilterFunc)loadSymbol(oidnLibrary, "oidnCommitFilter");
		oidn::ExecuteFilter = (oidn::ExecuteFilterFunc)loadSymbol(oidnLibrary, "oidnExecuteFilter");
		oidn::ReleaseFilter = (oidn::ReleaseFilterFunc)loadSymbol(oidnLibrary, "oidnReleaseFilter");
	}
};
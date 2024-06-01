#pragma once

NAMESPACE_BEGIN(Grid);

template<typename Vtype, typename obj>
class A2AViewBase {
public:
    std::vector<Vtype> _view;
    Vtype *_view_device;
    size_t _view_device_size;

public:
    A2AViewBase() = default;
    
    Vtype &operator[](size_t i) { return _view[i];}

    int size() {return _view.size();}

    void reserve(int size) { 
        _view.reserve(size); 
        _view_device_size = size*sizeof(Vtype);
        _view_device = (Vtype *)acceleratorAllocDevice(_view_device_size);
    }

    Vtype *getView() { return _view_device; }

    virtual void copyToDevice() {
        acceleratorCopyToDevice(_view.data(),_view_device,_view_device_size);
    }

    virtual void closeViews() {
        for(int p=0;p<this->_view.size();p++)   this->_view[p].ViewClose(); 

            _view.erase(this->_view.begin(),this->_view.end());

        if (this->_view_device_size > 0) {
            acceleratorFreeDevice(this->_view_device);
            this->_view_device_size = 0;
        }
    }

    virtual ~A2AViewBase() { this->closeViews();  }
};

template<typename obj>
class A2AFieldView: public A2AViewBase<LatticeView<obj>,obj> {
public:
    void openViews(const Lattice<obj> *fields, int size) {
        this->reserve(size);
        for (int i = 0; i < size; i++) {
            this->_view.push_back(fields[i].View(AcceleratorRead));
        }
        this->copyToDevice();
    }
};

template<typename obj,typename FImplParams>
class A2AStencilView: public A2AViewBase<CartesianStencilView<obj,obj,FImplParams>, obj>{
protected:
    std::vector<std::unique_ptr<CartesianStencil<obj,obj,FImplParams> > > _stencils;

    obj *_buffer_device;
    size_t _buffer_device_size;

    std::vector<Integer> _offset;
    Integer *_offset_device;
    size_t _offset_device_size;

public:
    obj *getBuffer() { return _buffer_device; }
    Integer *getOffset() { return _offset_device; }
    Integer & getOffset(int i) { return _offset[i]; }

    void addStencil(std::unique_ptr<CartesianStencil<obj,obj,FImplParams> > &stencil) {
        _stencils.push_back(std::move(stencil));
    }

    void openViews(const Lattice<obj> *fields, int size) {

        createCommBuffer(fields,size);

        this->reserve(_stencils.size());
        for (auto &stencil: _stencils) {
            this->_view.push_back(stencil->View(AcceleratorRead));
        }
        this->copyToDevice();
    }

    void createCommBuffer(const Lattice<obj> *fields, int nFields) {

        Vector<obj> buffer;

        GridBase *grid = fields[0].Grid();
        SimpleCompressor<obj> compressor;
        int comm_buf_size;
        obj *buf_p;

        int j = 0;
        bool multipleStencils = _stencils.size() > 1;
        if (multipleStencils) assert(_stencils.size() == nFields);

        for (int i = 0; i < nFields; i++) {
            const auto &field = fields[i];
            auto &stencil = _stencils[j];

            stencil->HaloExchange(field,compressor);

            comm_buf_size = stencil->_unified_buffer_size;
            _offset.push_back(buffer.size());

            buffer.resize(buffer.size()+comm_buf_size);
            buf_p = &buffer[_offset.back()];

            if (comm_buf_size > 0) {
                obj *comm_buf_p = stencil->CommBuf();
                accelerator_for (k, comm_buf_size, 1, {
                    buf_p[k] = comm_buf_p[k];
                });
            }

            if (multipleStencils) j++;
        }
        _offset_device_size = _offset.size()*sizeof(Integer);
        _offset_device = (Integer *)acceleratorAllocDevice(_offset_device_size);
        acceleratorCopyToDevice(_offset.data(),_offset_device,_offset_device_size);

        _buffer_device_size = buffer.size()*sizeof(obj);
        _buffer_device = (obj *)acceleratorAllocDevice(_buffer_device_size);
        acceleratorCopyToDevice(buffer.data(),_buffer_device,_buffer_device_size);

    }
    virtual void closeViews() {

        A2AViewBase<CartesianStencilView<obj,obj,FImplParams>, obj>::closeViews();
        
        _offset.resize(0);
        _stencils.resize(0);

        if (this->_buffer_device_size > 0) {
            acceleratorFreeDevice(_buffer_device);
            acceleratorFreeDevice(_offset_device);

            _offset_device_size = 0;
            _buffer_device_size = 0;
        }
    }

    ~A2AStencilView() { this->closeViews(); }
};

NAMESPACE_END(Grid);
